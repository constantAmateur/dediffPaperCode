###########################
# libraries & functions ##
##########################
library(dplyr)
library(plyr)
library(ggplot2)
library(data.table)
library(reshape2)
library(grid)
library(gridExtra)
source("/Users/mkalyva/GitHub/cellSignalAnalysis/cellSignalAnalysis_original.R") # script under cellSignalAnalysis tool available in github repository

# load tcga & hsc metadata
meta.tcga<-read.table('/Users/mkalyva/DATA/DECONVOLUTION/TCGA_Full_Metadata.tsv',header=T,sep="\t")
meta.tcga$UniqueSampleID<-gsub(x=meta.tcga$UniqueSampleID,pattern = "-",replacement = ".") # format metadata
# formating tcga and hsc metadata
meta.hsc<-read.table('/Users/mkalyva/DATA/DECONVOLUTION/blastoid_launch.tsv',header=F,sep="\t")
meta.hsc$V1 <- gsub(x = meta.hsc$V1,pattern = "tsv",replacement = "txt")

# define main signatures from early embryo & gastrulation
ee<-c("Hypoblast","Epiblast")
gast<-c("Endoderm","Mesoderm", "Non-Neural Ectoderm","Non.NeuralEctoderm")
int<-"Intercept"

# ------------------------------------------------------------------------------------------------------------------------------------
#  PART 1: load results from cellSignalAnalysis #
# ------------------------------------------------------------------------------------------------------------------------------------

#-- IMPORTANT NOTE ---#
# the data read in each time are the exposures from the cellSignalAnalysis tool for each reference used and the bulk-RNAseq data that was deconvoluted


# load early embryo
setwd("~/ee_fit/")
Data.in <- lapply(list.files(pattern = "\\.tsv$"),normaliseExposures)  
#combines all of the tables by row into one 
D<-do.call(cbind,as.data.frame(lapply(Data.in, `[[`, 1)))
rownames(D)<-row.names(Data.in[[1]]$exposures)
early_emb<-data.frame(Type="EarlyEmbryo",melt(D))

# load early embryo + gastrulation
setwd("~/TCGA_plot/ee_gast_fit/")
Data.in <- lapply(list.files(pattern = "\\.tsv$"),normaliseExposures)  
#combines all of the tables by row into one 
D<-do.call(cbind,as.data.frame(lapply(Data.in, `[[`, 1)))
rownames(D)<-row.names(Data.in[[1]]$exposures)
ee_gast<-data.frame(Type="EE&Gastrulation",melt(D))

# load tabula sapiens -adult-
setwd("~/TCGA_plot/ts_fit/")
Data.in <- lapply(list.files(pattern = "\\.tsv$"),normaliseExposures)  
#combines all of the tables by row into one 
D<-do.call(cbind,as.data.frame(lapply(Data.in, `[[`, 1)))
rownames(D)<-row.names(Data.in[[1]]$exposures)
ee_gast_ts<-data.frame(Type="EE&Gast&Adult",melt(D))


# load fetal(han) & adult(ts)
setwd("~TCGA_plot/Han_TS/")
Data.in <- lapply(list.files(pattern = "\\.tsv$"),normaliseExposures)  
#combines all of the tables by row into one 
D<-do.call(cbind,as.data.frame(lapply(Data.in, `[[`, 1)))
rownames(D)<-row.names(Data.in[[1]]$exposures)
ts_fetal <-data.frame(Type="EE&Gast& Fetal-Adult",melt(D))


######################
# merge data frames #
#####################

all_df <- rbind(early_emb,ee_gast,ee_gast_ts,ts_fetal) 
colnames(all_df) <- c("Type","celltype","sample","value") # add colnames

#########################
# necessary formatting ##
#########################

# when same cell type across tissues (eg. B cell) remove the end extension so that you do not consider them as different cell-types
all_df$celltype<- gsub(pattern = ".x",replacement = "",all_df$celltype, perl= F, fixed = T)
all_df$celltype<- gsub(pattern = ".y",replacement = "",all_df$celltype, perl = F, fixed = T)
all_df$celltype<- gsub(pattern = "_$",replacement = "",all_df$celltype, perl = T, fixed = T)

##############################################################################################################
## label cell-types under a bigger category: early embryo, gastrulation, fetal, adult or unexplained ########
#############################################################################################################
# ee, gast and int have been defined in first part of the script

# label cell types under early embryo, gastrulation first
all_df$celltype[which(all_df$Type=="EarlyEmbryo" & all_df$celltype %in% ee)]<-"EarlyEmbryo"
all_df$celltype[which(all_df$Type=="EE&Gastrulation" & all_df$celltype %in% ee)]<-"EarlyEmbryo"
all_df$celltype[which(all_df$Type=="EE&Gastrulation" & all_df$celltype %in% gast)]<-"Gastrulation"
all_df$celltype[which(all_df$Type=="EE&Gast&Adult" & all_df$celltype %in% ee)]<-"EarlyEmbryo"
all_df$celltype[which(all_df$Type=="EE&Gast&Adult" & all_df$celltype %in% gast)]<-"Gastrulation"
all_df$celltype[which(all_df$Type=="EE&Gast& Fetal-Adult" & all_df$celltype %in% ee)]<-"EarlyEmbryo"
all_df$celltype[which(all_df$Type=="EE&Gast& Fetal-Adult" & all_df$celltype %in% gast)]<-"Gastrulation"
#index cell-types that are only in adult reference (by removing all other refs)
idx_adult<-which(all_df$celltype != "EarlyEmbryo" & all_df$celltype != "Intercept" & all_df$celltype != "Gastrulation" & all_df$Type=="EE&Gast&Adult")
adult_types<-unique(all_df$celltype[idx_adult]) # save adult cell-types in variable
# label adult the cell-types from the adult reference
all_df$celltype[idx_adult] <- "Adult" 

# label adult cell-types from the joined adult-fetal reference 
all_df$celltype[which(all_df$Type=="EE&Gast& Fetal-Adult" & all_df$celltype %in% adult_types)]<-"Adult"

# label fetal cell-types from the joined adult-fetal reference 
idx_fetal<-which(all_df$celltype != "EarlyEmbryo" & all_df$celltype != "Intercept" & all_df$celltype != "Gastrulation" & all_df$celltype != "Adult" &  all_df$Type=="EE&Gast& Fetal-Adult")
fetal_types<-unique(all_df$celltype[idx_fetal])
all_df$celltype[idx_fetal] <- "Fetal"


#########################
# add tcga cancer type ##
##########################

# add a new layer on the table with the tissue origin
all_df$tcga<-"NULL"
tissues<-unique(meta.tcga$Tissue)


for(i in 1:length(tissues)){
  
  tis_now<-tissues[i]
  samples<-meta.tcga$UniqueSampleID[which(meta.tcga$Tissue==tis_now)]
  idx_now<-which(all_df$sample %in% samples)
  message(length(idx_now))
  all_df$tcga[idx_now]<-tis_now
  
}

all_df$tcga[grep("GSM5",all_df$sample)]<-"HSC" # add a tissue label for the control, the blastoids

# Add TCGA Abbreviations
meta.tcga$UniqueSampleID<-gsub("-",".",meta.tcga$UniqueSampleID)
idx.match<-match(all_df$sample,meta.tcga$UniqueSampleID)
all_df$main_types<-meta.tcga$main_types[idx.match]
all_df$main_types[which(all_df$tcga=="HSC")]<-"HSC"

# all_df table is the tcga data frame analysis

# ---- for plotting purposes only: factorise in the order you want the table results  --- #
# Type_f levels
all_df$Type_f = factor(all_df$Type, levels=c("EarlyEmbryo","EE&Gastrulation","EE&Gast&Adult","EE&Gast& Fetal-Adult"))
#celltype_f levels
c_levels<-c("Intercept","Fetal", "Adult","Gastrulation", "EarlyEmbryo")
all_df$celltype_f = factor(all_df$celltype, levels=c_levels)
#tc levels
tcga_levels<-unique(all_df$tcga)
tcga_levels<-tcga_levels[-9]
tcga_levels[20]<-"HSC" # put the control first
all_df$tcga_f = factor(all_df$tcga, levels=tcga_levels)
# --------------------------------------------------------------------------------------- #

# for each reference and sample aggregate 
aggr_all_df <- ddply( all_df, .(Type_f, celltype_f, tcga_f, sample), summarise, median_ref = sum(value))

# for each reference summarise mean and sd statistics
mean_aggr<-ddply( aggr_all_df, .(Type_f, celltype_f, tcga_f), summarise, sd_ref = sd(median_ref), median_ref = mean(median_ref))

# save analysis data
saveRDS(aggr_all_df, "~/tcga_aggr_all_df.RDS")
saveRDS(mean_aggr, "~/tcga_aggr_all_df_mean_sd.RDS")

#-------------------------
tcga_aggr_mean <- mean_aggr
tcga_aggr_mean$dataset <- "TCGA"
#-------------------------


# ------------------------------------------------------------------------------------------------------------------------------------
#  PART 2: load results from cellSignalAnalysis for pediatric cancers 
# ------------------------------------------------------------------------------------------------------------------------------------

#---------------------#
# READ in TARGET DATA #
#---------------------#
setwd("~/PED_Plot/exposures/")
# load early embryo
Data.in <- normaliseExposures("PED_EE_fitExposures.tsv")
#combines all of the tables by row into one 
D<-as.data.frame(Data.in$exposures)
D$cell_type<-rownames(D)
#rownames(D)<-row.names(Data.in[[1]]$exposures)
early_emb<-data.frame(Type="EarlyEmbryo",melt(D))

# load early embryo + gastrulation
Data.in <- normaliseExposures("PED_EE_GAST_fitExposures.tsv")
#combines all of the tables by row into one 
D<-as.data.frame(Data.in$exposures)
D$cell_type<-rownames(D)
ee_gast<-data.frame(Type="EE&Gastrulation",melt(D))

# load han - fetal-
Data.in <- normaliseExposures("PED_FETAL_fitExposures.tsv")
#combines all of the tables by row into one 
D<-as.data.frame(Data.in$exposures)
D$cell_type<-rownames(D)
ee_gast_fetal <-data.frame(Type="EE&Gast&Fetal",melt(D))


# load adult & fetal reference together
Data.in <- normaliseExposures("PED_FETAL_TS_fitExposures.tsv")
#combines all of the tables by row into one 
D<-as.data.frame(Data.in$exposures)
D$cell_type<-rownames(D)
ts_fetal <-data.frame(Type="EE&Gast& Fetal-Adult",melt(D))


# aggregate all references in one table
ped_df<-rbind(early_emb,ee_gast,ee_gast_fetal,ts_fetal)
colnames(ped_df)<-c("Type","celltype","sample","value") # add column names

# when same cell type across tissues (eg. B cell) remove the end extension so that you do not consider them as different cell-types
ped_df$celltype<- gsub(pattern = ".x",replacement = "",ped_df$celltype, perl= F, fixed = T)
ped_df$celltype<- gsub(pattern = ".y",replacement = "",ped_df$celltype, perl = F, fixed = T)

##############################################################################################################
## label cell-types under a bigger category: early embryo, gastrulation, fetal, adult or unexplained ########
#############################################################################################################
# ee, gast and int have been defined in first part of the script
# label cell types under early embryo, gastrulation first
ped_df$celltype[which(ped_df$Type=="EarlyEmbryo" & ped_df$celltype %in% ee)]<-"EarlyEmbryo"
ped_df$celltype[which(ped_df$Type=="EE&Gastrulation" & ped_df$celltype %in% ee)]<-"EarlyEmbryo"
ped_df$celltype[which(ped_df$Type=="EE&Gastrulation" & ped_df$celltype %in% gast)]<-"Gastrulation"
ped_df$celltype[which(ped_df$Type=="EE&Gast&Fetal" & ped_df$celltype %in% ee)]<-"EarlyEmbryo"
ped_df$celltype[which(ped_df$Type=="EE&Gast&Fetal" & ped_df$celltype %in% gast)]<-"Gastrulation"
ped_df$celltype[which(ped_df$Type=="EE&Gast& Fetal-Adult" & ped_df$celltype %in% ee)]<-"EarlyEmbryo"
ped_df$celltype[which(ped_df$Type=="EE&Gast& Fetal-Adult" & ped_df$celltype %in% gast)]<-"Gastrulation"

# index fetal cell-types from fetal reference
idx_fetal<-which(ped_df$celltype != "EarlyEmbryo" & ped_df$celltype != "Intercept" & ped_df$celltype != "Gastrulation" & ped_df$Type=="EE&Gast&Fetal")
fetal_cell_types<-unique(ped_df$celltype[idx_fetal])
ped_df$celltype[idx_fetal] <- "Fetal" # name fetal cell-types from the fetal reference

# name fetal cell-types from the joint fetal reference
ped_df$celltype[which(ped_df$Type=="EE&Gast& Fetal-Adult" & ped_df$celltype %in% fetal_cell_types)]<-"Fetal"

# name adult cell-types from the joint fetal reference
idx_adult<-which(ped_df$celltype != "EarlyEmbryo" & ped_df$celltype != "Intercept" & ped_df$celltype != "Gastrulation" & ped_df$celltype != "Fetal" &  ped_df$Type=="EE&Gast& Fetal-Adult")
adult_types<-unique(ped_df$celltype[idx_adult])
ped_df$celltype[idx_adult] <- "Adult"

##########################################################
### name the major types of pediatric data I have here ###
##########################################################

# remove sick kids dataset
idx_sick<-grep("Sick",ped_df$sample)
ped_df<-ped_df[-idx_sick,]

# rename cancer types from target
idx_WT<-grep("WT",ped_df$sample)
ped_df$target[idx_WT]<-"TARGET.WT"

idx_RT<-grep("RT",ped_df$sample)
ped_df$target[idx_RT]<-"TARGET.RT"

idx_AML<-grep("AML",ped_df$sample)
ped_df$target[idx_AML]<-"TARGET.AML"

idx_ALL<-grep("ALL",ped_df$sample)
ped_df$target[idx_ALL]<-"TARGET.ALL"

idx_NBL<-grep("NBL",ped_df$sample)
ped_df$target[idx_NBL]<-"TARGET.NBL"

# set levels 
ped_df$Type_f = factor(ped_df$Type, levels=c("EarlyEmbryo","EE&Gastrulation", "EE&Gast&Fetal","EE&Gast& Fetal-Adult"))
c_levels<-c("Intercept","Fetal", "Adult","Gastrulation", "EarlyEmbryo")
ped_df$celltype_f = factor(ped_df$celltype, levels=c_levels)

#tc levels
target_levels<-unique(ped_df$target)

aggr_ped_df <- ddply( ped_df, .(Type_f, celltype_f, target, sample), summarise, median_ref = sum(value))
# for each reference and sample aggregate 
mean_aggr_ped<-ddply( aggr_ped_df, .(Type_f, celltype_f, target), summarise, sd_ref = sd(median_ref), median_ref = mean(median_ref))

#-------------------------
TRG <- mean_aggr_ped
TRG$dataset <-"TARGET"
#-----------------------

#--------------------------#
# READ in St.Jude's DATA   #
#--------------------------#

# -- load StJudes -- #

meta.st.judes <- read.table("/Users/mkalyva/DATA/DECONVOLUTION/Pediatric/StJudesBasicMetadata.tsv", sep="\t", fill=T, header=T)
meta.st.judes.idx<-grep("RNA", meta.st.judes$sequencing_type)
keep<-c("ACT","BT","CBF","CPC","DSRCT","EPD","HGG","LGG","MB","MEL","NBL","RB","RHB","OS","ST","WLM")
meta.st.judes <- meta.st.judes[which(meta.st.judes$sj_diseases %in% keep),]

setwd("/Users/mkalyva/DATA/DECONVOLUTION/PED_Plot/exposures/")

# load early embryo
Data.in <- normaliseExposures("Judes_EE_fitExposures.tsv")
#combines all of the tables by row into one 
D<-as.data.frame(Data.in$exposures)
D$cell_type<-rownames(D)
#rownames(D)<-row.names(Data.in[[1]]$exposures)
early_emb<-data.frame(Type="EarlyEmbryo",melt(D))

# load load early embryo + gastrulation
#read in all .txt files but skip the first 8 rows
Data.in <- normaliseExposures("Judes_EE_GAST_fitExposures.tsv")
#combines all of the tables by row into one 
D<-as.data.frame(Data.in$exposures)
D$cell_type<-rownames(D)
ee_gast<-data.frame(Type="EE&Gastrulation",melt(D))

# load han et.al. - fetal -
Data.in <- normaliseExposures("Judes_Fetal_new_fitExposures.tsv")
#combines all of the tables by row into one 
D<-as.data.frame(Data.in$exposures)
D$cell_type<-rownames(D)
ee_gast_fetal <-data.frame(Type="EE&Gast&Fetal",melt(D))


# load  fetal& adult joint reference 
Data.in <- normaliseExposures("Judes_TS_Fetal_new_fitExposures.tsv")
#combines all of the tables by row into one 
D<-as.data.frame(Data.in$exposures)
D$cell_type<-rownames(D)
ts_fetal <-data.frame(Type="EE&Gast& Fetal-Adult",melt(D))

ped_df<-rbind(early_emb,ee_gast,ee_gast_fetal,ts_fetal)
colnames(ped_df)<-c("Type","celltype","sample","value")
ped_df$celltype<- gsub(pattern = ".x",replacement = "",ped_df$celltype, perl= F, fixed = T)
ped_df$celltype<- gsub(pattern = ".y",replacement = "",ped_df$celltype, perl = F, fixed = T)

##############################################################################################################
## label cell-types under a bigger category: early embryo, gastrulation, fetal, adult or unexplained ########
#############################################################################################################
# ee, gast and int have been defined in first part of the script
# label cell types under early embryo, gastrulation first
ped_df$celltype[which(ped_df$Type=="EarlyEmbryo" & ped_df$celltype %in% ee)]<-"EarlyEmbryo"
ped_df$celltype[which(ped_df$Type=="EE&Gastrulation" & ped_df$celltype %in% ee)]<-"EarlyEmbryo"
ped_df$celltype[which(ped_df$Type=="EE&Gastrulation" & ped_df$celltype %in% gast)]<-"Gastrulation"
ped_df$celltype[which(ped_df$Type=="EE&Gast&Fetal" & ped_df$celltype %in% ee)]<-"EarlyEmbryo"
ped_df$celltype[which(ped_df$Type=="EE&Gast&Fetal" & ped_df$celltype %in% gast)]<-"Gastrulation"
ped_df$celltype[which(ped_df$Type=="EE&Gast& Fetal-Adult" & ped_df$celltype %in% ee)]<-"EarlyEmbryo"
ped_df$celltype[which(ped_df$Type=="EE&Gast& Fetal-Adult" & ped_df$celltype %in% gast)]<-"Gastrulation"

# index fetal from the fetal reference
idx_fetal<-which(ped_df$celltype != "EarlyEmbryo" & ped_df$celltype != "Intercept" & ped_df$celltype != "Gastrulation" & ped_df$Type=="EE&Gast&Fetal")
fetal_cell_types<-unique(ped_df$celltype[idx_fetal])
ped_df$celltype[idx_fetal] <- "Fetal" 
# name fetal from the fetal & adult joint reference
ped_df$celltype[which(ped_df$Type=="EE&Gast& Fetal-Adult" & ped_df$celltype %in% fetal_cell_types)]<-"Fetal"
# name adult from the fetal & adult joint reference
idx_adult<-which(ped_df$celltype != "EarlyEmbryo" & ped_df$celltype != "Intercept" & ped_df$celltype != "Gastrulation" & ped_df$celltype != "Fetal" &  ped_df$Type=="EE&Gast& Fetal-Adult")
adult_types<-unique(ped_df$celltype[idx_adult])
ped_df$celltype[idx_adult] <- "Adult"

samples<-unique(meta.st.judes$sample_name)
# remove st.jude's cancer types not relevant to adult or too low in samples
ped_df<-ped_df[which(ped_df$sample %in% samples),]
# set levels for facet grid
ped_df$Type_f = factor(ped_df$Type, levels=c("EarlyEmbryo","EE&Gastrulation", "EE&Gast&Fetal","EE&Gast& Fetal-Adult"))
c_levels<-c("Intercept","Fetal", "Adult","Gastrulation", "EarlyEmbryo")
ped_df$celltype_f = factor(ped_df$celltype, levels=c_levels)
#order sample levels
ped_df$sample_order<-"NULL"
idx<-match(ped_df$sample, meta.st.judes$sample_name)
ped_df$sample_order<-meta.st.judes$sj_diseases[idx]

ped_df <- ped_df[-grep(pattern = ".1$",perl=T, x = ped_df$sample),]
ped_df <- ped_df[-grep(pattern = ".2$",perl=T, x = ped_df$sample),]
ped_df <- ped_df[-grep(pattern = ".3$",perl=T, x = ped_df$sample),]
ped_df <- ped_df[-grep(pattern = ".4$",perl=T, x = ped_df$sample),]
#ped_df <- ped_df[-grep(pattern = ".5$",perl=T, x = ped_df$sample),]

# aggregate tissue specific cell-signals
# for each reference and sample aggregate 
aggr_ped_df <- ddply( ped_df, .(Type_f, celltype_f, sample_order, sample), summarise, median_ref = sum(value))
# for each reference and sample aggregate 
mean_aggr_ped<-ddply( aggr_ped_df, .(Type_f, celltype_f, sample_order), summarise, sd_ref = sd(median_ref), median_ref = mean(median_ref))
#--------------------------
Judes_mean <- mean_aggr_ped
Judes_mean$dataset <- "St.Jude"
#--------------------------


# ------------------------------------------------------------------------------------------------------------------------------------
#  PART 3: Merge results from adult & pediatric cancers and plot together (Figure 1)
# ------------------------------------------------------------------------------------------------------------------------------------

# -------- read and merge TCGA, TARGET & StJude's results all together ------------ #
colnames(TRG)<-colnames(Judes_mean)
colnames(tcga_aggr_mean)[3]<-"sample_order"
df <- rbind(tcga_aggr_mean,TRG,Judes_mean)

# # order samples by ee embryo signal -- for plotting purposes only #
# df2<-mean_aggr[which(mean_aggr$celltype_f=="EarlyEmbryo" & mean_aggr$Type_f=="EarlyEmbryo"),]
# df_sorted <- df2 %>% 
#   arrange(median_ref)
# sample_order<-unique(df_sorted$tcga_f)
# mean_aggr$tcga_order<-factor(mean_aggr$tcga_f, levels=sample_order)

# add levels to organise order in plotting
df$dataset <- as.character(df$dataset)
df$sample_order <- as.character(df$sample_order)
df$sample_order[which(df$sample_order=='HSC')] <- "Blastoids"
df$dataset[which(df$sample_order=='Blastoids')] <- "Blastoids"
df$dataset_f <- factor(df$dataset, levels=c("Blastoids","TCGA","TARGET","St.Jude"))


df2<-df[which(df$celltype_f=="EarlyEmbryo" & df$Type_f=="EarlyEmbryo" & df$dataset != "St.Jude"),]
df_sorted <- df2 %>% 
  arrange(median_ref)
ee_sample_order<-unique(df_sorted$sample_order) # orders based on ee signal

df2<-df[which(df$celltype_f=="EarlyEmbryo" & df$Type_f=="EarlyEmbryo" & df$dataset == "St.Jude"),]
df_sorted <- df2 %>% 
  arrange(median_ref)
judes_order<-unique(df_sorted$sample_order) # order of st.jude

ee_sample_order <- ee_sample_order[-(grep("TARGET",ee_sample_order))] # for target do not order based on ee signal alone
ee_sample_order <- as.character(ee_sample_order[-(grep("Blastoids",ee_sample_order))])
all_sample_order <- c("Blastoids", ee_sample_order,"TARGET.NBL","TARGET.WT","TARGET.RT", judes_order) # final order of cancers for the plot

# remove leucemias from TARGET
leu<-c("TARGET.ALL","TARGET.AML")
df <- df[-which(df$sample_order %in% leu),]

# ----- plot ------ #

df$ee_sample_order <- factor(df$sample_order, levels=all_sample_order)

df <- df[-which(df$Type_f=="EE&Gast&Adult"),]
df <- df[-which(df$Type_f=="EE&Gast&Fetal"),]
#df <- df[-which( is.na(df$ee_sample_order)),] 
df <- df[-which(df$sample_order=="CBF"),] # remove St Judes CBF

# Set colours for barplot
myColors <- c("#767676FF","#8888FFFF","#CC79A7FF","#8F3931FF","#4DBBD5FF")
colScale <- scale_fill_manual("Cell type",values = myColors,drop=FALSE)

# plot pancancer barplot

pdf(file = "/Users/mkalyva/DATA/DECONVOLUTION/Figure_1_by_Tissue_Pediatric_and_TCGA.pdf",width = 32, height = 18 ,onefile = T)
ggplot(data = df,
       mapping = aes(x = ee_sample_order, y = median_ref, fill= celltype_f)) +
  colScale+ ylab("Cellular Signal") + xlab("Cancer Type")+
  geom_bar(stat="identity", position = "stack")+
  theme_bw() + scale_y_continuous(position = "left") + ylab("\nCellular Signal")+
  facet_grid(dataset_f ~ Type_f, scales = "free", space = "free")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        text = element_text(size=16),axis.text.x=element_text(angle=90, hjust=1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 10, l = 10)),
        strip.text.y = element_text(angle=0), strip.background = element_rect(colour="black",fill='white',
                                                                              size=1.5, linetype="solid")) +
  coord_flip()+
  guides(fill=guide_legend(title="Cell Type")) # + geom_errorbar(limits,width=.2,position=dodge, size=0.1)  # Width of the error bars
dev.off()




