#' Create (some of) the figure panels relating to the analysis of the liver data.  Based on Gerda's code and data...
setwd('~/trueHome/Projects/DeDiff')

#############
# Libraries #
#############

library(SoupX)
library(Seurat)
library(ggplot2)
library(reshape2)
library(stringr)
library(ComplexHeatmap)
library(Matrix)
#import('Code/cellTypist.R',as='ct')
import('Code/logisticRegressionCellTypist.R',as='ct')
import('Code/cellSignalAnalysis.R',as='csa')


##########
# Params #
##########

nParallel=48
plotDir = 'Results/GutFigures/'
s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes
srcGut = 'Results/rasa_rel_srat.rds'
srcGutGenes = 'Results/gutGenes.csv'
srcAdr = 'Results/adr_all.rds'
srcAdrGenes = 'Results/adrGenes.tsv.gz'
srcGutWithEE = 'Results/gut_gast_ee_ref.rds'
srcCRC = 'Data/LeeCRC/'
srcGutCS = 'Results/cellSignalOutput/Gerdas/gutFull'
srcBoardman = 'Data/boardmanData/'
fitModelPaths = c(gutGastEE = 'Results/gutGastEE_CT.pkl')

#############
# Functions #
#############

process_seurat =function(srat){
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat)
  srat = RunPCA(srat, npcs = 50)
  srat = FindNeighbors(srat, dims=1:50)
  srat = FindClusters(srat, resolution = 1)
  srat = RunUMAP(srat, dims=1:50, min.dist = 0.5, n.neighbors = 50)
  return(srat)
}

#######################
# Cell signal summary #
#######################

plotType='boxplotsByExposure'
#Load the main fits and some metadata
gut_cellsig=csa$normaliseExposures(file.path(srcGutCS,'combinedCortex_OutRun_fitExposures.tsv'))
paths=read.table(file.path(srcGutCS,'all_gut_paths.txt'), header = F)
paths$dataset=c(rep("TCGA", 715), rep("Boardman", 79), rep("fetal bulk", 8), rep("gtex_brain", 370))
paths$sample=basename(paths$V1)
#Boardman metadata 
boardman_meta=read.table(file.path(srcBoardman,'boardman_meta.txt'), sep = ",", header = T)
boardman_meta_new=data.frame(sample=colnames(gut_cellsig$exposures)[716:794])
#which(boardman_meta_new$sample%in%boardman_meta$sample)
boardman_meta_new$sample_type="x"
boardman_meta_new$sample_type=boardman_meta$sample_type[match(boardman_meta_new$sample, boardman_meta$sample)]
boardman_meta_new$sample_type[which(is.na(boardman_meta_new$sample_type))]="unknown"
boardman_meta_new$sample_type_wide=boardman_meta_new$sample_type
boardman_meta_new$sample_type_wide[grep("VILLOUS", boardman_meta_new$sample_type_wide)]= "VILLOUS ADENOMA"
#More metadata
mDat_col=readRDS(file.path(srcGutCS,'mDat_gut.rds'))
#TCGA
tcga_meta=data.frame(sample=colnames(gut_cellsig$exposures)[1:715],
                     sample_type=c(as.character(mDat_col$gdc_cases.samples.sample_type)),
                     stage=c(mDat_col$xml_stage_event_pathologic_stage))
paths$sample_type="x"
paths$sample_type[1:715]=tcga_meta$sample_type
paths$sample_type[716:794]=boardman_meta_new$sample_type_wide
paths$sample_type[795:802]="fetal_bulk"
paths$sample_type[803:nrow(paths)]="gtex_brain"
#Load full thing just to plot as meta-data environment is right for plot here
pdf(file.path(plotDir,'cellSignalGutFull.pdf'),width=9,height=7)
paths$labs = c(CANCER='Colorectal cancer',
               'NORMAL EPITH' = 'Normal adult gut',
               'Solid Tissue Normal' = 'Normal adult gut',
               'VILLOUS ADENOMA' = 'Villous adenoma',
               'fetal_bulk' = 'Normal fetal gut',
               'gtex_brain' = 'GTEX brain',
               'Primary Tumor' = 'Colorectal cancer')[paths$sample_type]
t1=csa$normaliseExposures(file.path(srcGutCS,'everything_OutRun_fitExposures.tsv'))
w = which(!is.na(paths$labs))
t2 = lapply(t1,function(e) e[,w])
t2$exposures = t2$exposures[c('Endoderm','Mesoderm','Non-Neural Ectoderm','Epiblast','Hypoblast',grep('Fetal_',rownames(t2$exposures),value=TRUE),grep('Adult_',rownames(t2$exposures),value=TRUE),'Cortex'),]
hm = csa$plotExposures(t2, 
                       column_split=paths$labs[w], 
                       show_column_names=FALSE, 
                       show_column_dend=FALSE,
                       column_title_rot=90, 
                       cluster_column_slices=FALSE,
                       cluster_rows=FALSE,
                       height = nrow(gut_cellsig$exposures)*unit(3, "mm"), column_gap=unit(3, "mm"))
draw(hm)
dev.off()
master_df=as.data.frame(t(gut_cellsig$exposures))
master_df$sample=rownames(master_df)
master_df$total_adult_signal=rowSums(master_df[,grep("Adult", colnames(master_df))])
master_df$total_fetal_signal=rowSums(master_df[,grep("Fetal", colnames(master_df))])
master_df$Adult_stem_ta=(master_df$`Adult_Stem cells`+master_df$Adult_TA)
master_rel_cols=master_df[,c("sample", "Intercept", "total_fetal_signal", "total_adult_signal", "Adult_stem_ta", "Cortex")]
master_rel_cols$sample_type=paths$sample_type
master_rel_cols$dataset=paths$dataset
#table(master_rel_cols$sample_type)
master_rel_cols$sample_type[which(master_rel_cols$sample_type=="CANCER")]="Primary Tumor"
master_rel_cols$sample_type[which(master_rel_cols$sample_type=="NORMAL EPITH")]="Solid Tissue Normal"
master_rel_cols$sample_type[which(master_rel_cols$sample_type=="VILLOUS ADENOMA")]="Villous Adenoma"
master_rel_cols=master_rel_cols[master_rel_cols$sample_type%in%c("fetal_bulk","gtex_brain",
                                                 "Primary Tumor", "Solid Tissue Normal", "Villous Adenoma" ),]
master_melted=melt(master_rel_cols, id.vars = c("sample", "sample_type", "dataset"))
master_melted$sample_type=factor(master_melted$sample_type, levels = c("gtex_brain","fetal_bulk", "Solid Tissue Normal",
                                                                       "Villous Adenoma","Primary Tumor"))
master_melted$dataset[master_melted$dataset=="fetal bulk"]="Kraiczy (Fetal bulk)"
master_melted$dataset[master_melted$dataset=="gtex_brain"]="GTEx (Brain)"
master_melted$dataset=factor(master_melted$dataset, levels=c("Kraiczy (Fetal bulk)","Boardman","TCGA", "GTEx (Brain)"))

####################
# Complete analysis
#The above analysis is preliminary, but good for defining metadata
gut_meta=read.table(file.path(srcGutCS,'gut_bulk_sample_metadata.txt'), sep = "\t", header = T)
rownames(gut_meta)=gut_meta$sample_name_no_dashes
gut_cellsig=csa$normaliseExposures(file.path(srcGutCS,'everything_OutRun_fitExposures.tsv'))
#mDat
mDat_gut=readRDS(file.path(srcGutCS,'mDat_gut.rds'))
gut_cellsig_df=as.data.frame(t(gut_cellsig$exposures))
#gut_cellsig_df$best_lab=colnames(gut_cellsig_df)[apply(gut_cellsig_df,1,which.max)]
#Replicate the meta-data pull from above
master_df=as.data.frame(t(gut_cellsig$exposures))
master_df$sample=rownames(master_df)
master_df$early_embryo = master_df$Epiblast + master_df$Hypoblast
master_df$gastrulation = master_df$Endoderm + master_df$Mesoderm + master_df$`Non-Neural Ectoderm`
master_df$fetal = rowSums(master_df[,grep('^Fetal_',colnames(master_df))])
master_df$adult = rowSums(master_df[,grep('^Adult_',colnames(master_df))])
master_df$adult_stemEpithelium = master_df$`Adult_Stem cells`+master_df$Adult_TA
master_df$adult_matureEpithelium = rowSums(master_df[,c('Adult_BEST4+ epithelial','Adult_Colonocyte')])
master_df$adult_specialisedEpithelium = rowSums(master_df[,c('Adult_Tuft','Adult_Goblet cell','Adult_Paneth','Adult_Microfold cell','Adult_BEST2+ Goblet cell')])
master_df$adult_other = rowSums(master_df[,grep('^Adult_',colnames(master_df))])-rowSums(master_df[,grep('adult_',colnames(master_df))])
master_df$adult_otherEpithelium = rowSums(master_df[,c('adult_matureEpithelium','adult_specialisedEpithelium')])
master_df$negative_control = master_df$Cortex
#master_df$total_adult_signal=rowSums(master_df[,grep("Adult", colnames(master_df))])
#master_df$total_fetal_signal=rowSums(master_df[,grep("Fetal", colnames(master_df))])
#master_df$Adult_stem_ta=(master_df$`Adult_Stem cells`+master_df$Adult_TA)
master_df$sample_type = paths$sample_type
master_df$dataset = paths$dataset
master_df$sample_type[which(master_df$sample_type=="CANCER")]="Primary Tumor"
master_df$sample_type[which(master_df$sample_type=="NORMAL EPITH")]="Solid Tissue Normal"
master_df$sample_type[which(master_df$sample_type=="VILLOUS ADENOMA")]="Villous Adenoma"
master_df=master_df[master_df$sample_type%in%c("fetal_bulk","gtex_brain",
                                                 "Primary Tumor", "Solid Tissue Normal", "Villous Adenoma" ),]
if(TRUE || plotType=='boxplotsByExposure'){
  dd = master_df[,c('early_embryo','gastrulation','fetal','adult_stemEpithelium','adult_otherEpithelium','Intercept','negative_control','sample_type')]
}else{
  #Really the only difference between this and the above is explicitly adding in the extra non-epithelial adult population to ensure things sum to 1
  #dd = master_df[,c('early_embryo','gastrulation','fetal','adult_stemEpithelium','adult_otherEpithelium','Intercept','negative_control','sample_type')]
  dd = master_df[,c('early_embryo','gastrulation','fetal','adult_stemEpithelium','adult_matureEpithelium','adult_specialisedEpithelium','adult_other','Intercept','negative_control','sample_type')]
}
tmp = table(dd$sample_type)
#Define ordering
tmp = tmp[c('gtex_brain','fetal_bulk','Solid Tissue Normal','Villous Adenoma','Primary Tumor')] 
#To generate the adult colours, use colorRampPalette around this colour '#CC79A7'
colfunc <- colorRampPalette(c('#E0C6D5', "#7F1153"))
colMap = c('Intercept'='lightgrey',
           'early_embryo'='#4DBBD5',
           'gastrulation'='#8F3931',
           'fetal'='#8888FF',
           'adult'='#CC79A7',
           'adult_stemEpithelium'="#E0C6D5",
           'adult_matureEpithelium'="#BF89A9",
           'adult_specialisedEpithelium'="#9F4D7E",
           'adult_other'="#7F1153",
           'negative_control'='black')
pdf(file.path(plotDir,'cellSignalGut.pdf'),width=7,height=3.5)
if(plotType=='boxplotsByExposure'){
  layout(matrix(seq(ncol(dd)-1),nrow=1),widths=c(1.3,rep(1,ncol(dd)-2)))
  for(tgt in colnames(dd[-ncol(dd)])){
    ii = match(tgt,colnames(dd))
    if(ii==1){
      par(mar=c(5,3,1,0.5))
    }else{
      par(mar=c(5,1,1,0.5))
    }
    s = split(dd[,tgt],dd$sample_type)[names(tmp)]
    plot(0,
         xlab='',
         ylab='exposure',
         main=tgt,
         xaxt='n',
         yaxt='n',
         xlim=c(0,length(s)),
         ylim=c(0,1),
         type='n')
    if(ii==1)
      axis(2,at=c(0,0.5,1),label=c(0,0.5,1))
    for(i in seq_along(s)){
      x = s[[i]]
      #Calibrate transparency so dots don't overwhelm the boxplot.  That is alpha*nPonits = pointsToShow
      alpha = min(1,30/length(x))
      alpha = toupper(as.hexmode(round(255*alpha)))
      alpha = paste0('#',paste(rep('0',6+(nchar(alpha)==1)),collapse=''),alpha)
      points(i-0.5+rnorm(length(x),sd=0.1),
             y=x,
             pch=19,
             col=alpha,
             cex=0.1)
      lines(x=c(i-0.5-0.3,i-0.5+0.3),
           y=rep(median(x),2),
           lwd=2)
      lines(x=rep(i-0.5,2),
           y=quantile(x,c(0.25,0.75)),
           lwd=2)
    }
    axis(1,at=seq(0.5,length(s)-0.5),label=names(s),las=2)
  }
}else{
  layout(matrix(seq_along(tmp),nrow=1))
  par(mar=c(5,3,1,1))
  for(dt in names(tmp)){
    x = dd[dd$sample_type==dt,-ncol(dd)]
    if(plotType=='stackedBars'){
      #Find whichever is the maximum
      m = which.max(colMeans(x))
      #Sort by the sum of the ones up to just before the max 
      x = x[order(rowSums(x[,1:max(1,(m-1)),drop=FALSE])),]
      x = x[hclust(dist(x))$order,]
      barplot(t(x),
              border=NA,
              space=0,
              yaxt='n',
              xaxt='n',
              col=colMap[colnames(x)])
      title(xlab=dt,line=0)
    }else if(plotType=='boxplots'){
      #Place each contribution at the same position on the x-axis then jitter dots around it horizontally
      plot(0,
           xlab='',
           ylab='exposure',
           main=dt,
           xaxt='n',
           xlim=c(0,ncol(x)),
           ylim=c(0,1),
           type='n')
      for(i in seq(ncol(x))){
        points(i-0.5+rnorm(nrow(x),sd=0.1),
               y=x[,i],
               pch=19,
               col='#00000034',
               cex=0.1)
        lines(x=c(i-0.5-0.3,i-0.5+0.3),
             y=rep(median(x[,i]),2),
             lwd=2)
      }
      axis(1,at=seq(0.5,ncol(x)-0.5),label=colnames(x),las=2)
    }
  }
}
dev.off()


####################
# Create reference #
####################

gut_srat=readRDS(srcGut)
# create annotation where epithelial cells have narrow annotation, but everything else has broad annotation
gut_srat@meta.data$annotation_narrow=gut_srat@meta.data$Integrated_05
gut_srat@meta.data$annotation_narrow[which(gut_srat@meta.data$annotation_narrow%in%c("D cells (SST+)","EC cells (NPW+)",
                                                                                     "EC cells (TAC1+)",
                                                                                     "EECs","I cells (CCK+)",
                                                                                     "L cells (PYY+)","M/X cells (MLN/GHRL+)",
                                                                                     "N cells (NTS+)",
                                                                                     "Progenitor (NEUROG3+)"))]="Enteroendocrine"
gut_srat@meta.data$annotation_wide=gut_srat@meta.data$category
gut_srat@meta.data$fetal_or_adult=gut_srat@meta.data$Diagnosis
# distinguish fetal and adult cells
gut_srat@meta.data$fetal_or_adult[which(gut_srat@meta.data$fetal_or_adult=="fetal")]="Fetal_"
gut_srat@meta.data$fetal_or_adult[which(gut_srat@meta.data$fetal_or_adult=="Healthy adult")]="Adult_"
gut_srat@meta.data$annotation_wide=paste0(gut_srat@meta.data$fetal_or_adult, gut_srat@meta.data$annotation_wide)
gut_srat@meta.data$annotation_narrow=paste0(gut_srat@meta.data$fetal_or_adult,gut_srat@meta.data$annotation_narrow)
gut_srat@meta.data$final_annot=gut_srat@meta.data$annotation_wide
gut_srat@meta.data$final_annot[which(gut_srat@meta.data$final_annot%in%c("Adult_Epithelial", "Fetal_Epithelial"))]=gut_srat@meta.data$annotation_narrow[which(gut_srat@meta.data$final_annot%in%c("Adult_Epithelial", "Fetal_Epithelial"))]
#remove clusters with <10 cells 
#table(gut_srat$final_annot)
gut_srat=subset(gut_srat, final_annot%in%c("Fetal_CLDN10+ cells","Fetal_Microfold cell", "Fetal_Paneth", "Fetal_Tuft"), invert=T)
#subset cortex from adrenal gland dataset
adr_all=readRDS(srcAdr)
cortex=subset(adr_all, Annotation%in%c("Cortex"))
cortex@meta.data$final_annot=cortex@meta.data$Annotation
#relabel the gene names to match
gut_genes=read.csv(srcGutGenes)
adr_genes=read.table(srcAdrGenes, sep = "\t", header = F)
gut_mtx=gut_srat@assays$RNA@counts
cortex_mtx=cortex@assays$RNA@counts
rownames(gut_mtx)=gut_genes$gene_ids
rownames(cortex_mtx)=adr_genes$V1
inter=intersect(rownames(gut_mtx), rownames(cortex_mtx))
rownames(gut_genes)=gut_genes$gene_ids
gut_mtx=gut_mtx[inter, ]
cortex_mtx=cortex_mtx[inter, ]
inter_genes=gut_genes[inter, ]
rownames(adr_genes)=adr_genes$V1
inter_adr=adr_genes[inter, ]
rownames(gut_mtx)=inter_adr$V2
rownames(cortex_mtx)=inter_adr$V2
new_g_srat=CreateSeuratObject(gut_mtx, meta.data = gut_srat@meta.data)
new_g_cortex=CreateSeuratObject(cortex_mtx, meta.data = cortex@meta.data)
new_g_merged=merge(new_g_srat, new_g_cortex)


############
# Load CRC #
############

#lee colorec data
crc_lee=read.table(file.path(srcCRC,'GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt'), header = T, sep = "\t", row.names = 1)
crc_lee_annot=read.table(file.path(srcCRC,'GSE132465_GEO_processed_CRC_10X_cell_annotation.txt'), header = T, sep = "\t")
crc_lee=as.matrix(crc_lee)
crc_lee=Matrix(crc_lee, sparse=T)
crc_lee_annot$Index=str_replace(crc_lee_annot$Index, pattern = "-", replacement = ".")
crc_lee_annot_tum=crc_lee_annot[which(crc_lee_annot$Class=="Tumor"),]
#dim(crc_lee_annot)
crc_lee_tum=crc_lee[,crc_lee_annot_tum$Index]
rownames(crc_lee_annot_tum)=crc_lee_annot_tum$Index
crc_lee_srat=CreateSeuratObject(crc_lee_tum, meta.data = crc_lee_annot_tum)
crc_lee_srat=process_seurat(crc_lee_srat)
crc_lee_srat@meta.data$Cell_subtype
crc_lee_srat@meta.data$wide_ann_pt=crc_lee_srat@meta.data$Cell_type
crc_lee_srat@meta.data$wide_ann_pt[which(crc_lee_srat@meta.data$wide_ann_pt=="Epithelial cells")]=paste0("Epithelial_cells_", crc_lee_srat@meta.data$Patient[which(crc_lee_srat@meta.data$wide_ann_pt=="Epithelial cells")])
crc_lee_srat@meta.data$wide_ann_pt[which(crc_lee_srat@meta.data$wide_ann_pt=="Epithelial cells")]=paste0(crc_lee_srat@meta.data$Cell_subtype[which(crc_lee_srat@meta.data$wide_ann_pt=="Epithelial cells")], crc_lee_srat@meta.data$Patient[which(crc_lee_srat@meta.data$wide_ann_pt=="Epithelial cells")])
#saveRDS(crc_lee_srat,'/lustre/scratch117/casm/team274/gk14/Dediff/GSE132465/crc_lee_srat.rds')
#crc_lee_srat=readRDS('/lustre/scratch117/casm/team274/gk14/Dediff/GSE132465/crc_lee_srat.rds')
#Simplify even furhter
crc_lee_srat@meta.data$annotPlot = crc_lee_srat@meta.data$wide_ann_pt
crc_lee_srat@meta.data$annotPlot[which(crc_lee_srat@meta.data$annotPlot %in% c('T cells','Mast cells','B cells','Myeloids'))]='Leukocytes'
crc_lee_srat@meta.data$annotPlot[as.character(crc_lee_srat@active.ident)=='22']='Endothelium'
colfunc <- colorRampPalette(c("#084594", "#C6DBEF"))
pdf(file.path(plotDir,"crc_umap.pdf"), height = 3, width=3,useDingbats = F)
gg=DimPlot(crc_lee_srat, group.by = "annotPlot", label=T,
        #cols = c( "#771122",colfunc(23),"#771122","#DDCC77","#44AA99","#7E481C"),
        cols = c('#A5636D', colfunc(23),"#DDCC77","#44AA99"),
        pt.size = 0.1) +theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
                              axis.ticks = element_blank(), legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1") +ggtitle("")
plot(gg)
dev.off()

#######################
# Logistic regression #
#######################

#############
# Fit models
#Need to load this thing, which will be the full reference including the EE stuff
gut_gastr_ee = readRDS(srcGutWithEE)
#tum_inter=intersect(rownames(new_g_merged), rownames(crc_lee_srat))
tum_ee_inter=intersect(rownames(gut_gastr_ee), rownames(crc_lee_srat))
#Now train the model
tgtFile = fitModelPaths['gutGastEE']
if(!file.exists(tgtFile)){
  ct$trainCelltypistModel(gut_gastr_ee@assays$RNA@counts[tum_ee_inter,],
                          outPath = tgtFile,
                          n_jobs=nParallel,
                          labels = gut_gastr_ee@meta.data$annot)
}
pr_fullref = ct$runCelltypist(crc_lee_srat@assays$RNA@counts[tum_ee_inter,],
                              model=tgtFile)

##########################
# Define column groupings
gut_colsplit_df=data.frame(celltype=colnames(pr_fullref$logitMat))
gut_colsplit_df$wide_annot=gut_colsplit_df$celltype
gut_colsplit_df$wide_annot[which(gut_colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Hypoblast", "Epiblast")))]="Early embryo"
gut_colsplit_df$wide_annot[which(gut_colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Mesoderm",
                                                                 "Endoderm",
                                                                 "Non-Neural Ectoderm")))]="Gastrulation"
gut_colsplit_df$wide_annot[which(gut_colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Fetal_BEST4+ epithelial",
                                                                 "Fetal_Colonocyte","Fetal_Distal progenitor", 
                                                                 "Fetal_Enterocyte",
                                                                 "Fetal_Goblet cell", 
                                                                 "Fetal_Proximal progenitor",
                                                                 "Fetal_Stem cells" ,
                                                                 "Fetal_TA")))]="Fetal_epithelial"
gut_colsplit_df$wide_annot[which(gut_colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Adult_Paneth",
                                                                 "Adult_BEST2+ Goblet cell",
                                                                 "Adult_BEST4+ epithelial",
                                                                 "Adult_Colonocyte", 
                                                                 "Adult_Goblet cell",
                                                                 "Adult_Microfold cell",
                                                                 "Adult_Stem cells",
                                                                 "Adult_Tuft",
                                                                 "Adult_TA")))]="Adult_epithelial"
gut_colsplit_df$wide_annot[which(gut_colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Fetal_Mesenchymal")))]='Fetal_Fibroblasts'
gut_colsplit_df$wide_annot[which(gut_colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Adult_Mesenchymal")))]='Adult_Fibroblasts'
gut_colsplit_df$wide_annot[which(gut_colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Adult_B cells",
                                                                 "Adult_Myeloid",
                                                                 "Adult_Plasma cells",
                                                                 "Adult_T cells")))]="Adult_leukocytes"
gut_colsplit_df$wide_annot[which(gut_colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Fetal_B cells" ,
                                                                 "Fetal_T cells", "Fetal_Myeloid" )))]="Fetal_leukocytes"
gut_colsplit_df$wide_annot[which(gut_colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Fetal_Red blood cells" )))]="Other blood cells"
gut_colsplit_df$wide_annot[which(gut_colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Adult_Enteroendocrine",
                                                                 "Adult_Mesenchymal",
                                                                 "Adult_Neuronal" )))]="Adult_Fibroblasts"
gut_colsplit_df$wide_annot[which(gut_colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Fetal_Enteroendocrine",
                                                                 "Fetal_Neuronal",
                                                                 "Fetal_Mesenchymal")))]="Fetal_Fibroblasts"
gut_colsplit_df$wide_annot=factor(gut_colsplit_df$wide_annot, levels=c("Early embryo" , 
                                                                       "Gastrulation" ,
                                                                       "Adult_epithelial" ,
                                                                       "Fetal_epithelial",
                                                                       'Adult_Endothelial',
                                                                       'Fetal_Endothelial',
                                                                       "Adult_leukocytes",
                                                                       "Fetal_leukocytes",
                                                                       'Other blood cells',
                                                                       "Adult_Fibroblasts",
                                                                       "Fetal_Fibroblasts",
                                                                       "Cortex"))
gut_colsplit_df$celltype=factor(gut_colsplit_df$celltype, levels=gsub('[ -]|\\+','.',c(
  "Hypoblast", "Epiblast", "Mesoderm",
  "Endoderm",
  "Non-Neural Ectoderm",
  "Adult_Stem cells",
  "Adult_TA",
  "Adult_BEST2+ Goblet cell",
  "Adult_Paneth",
  "Adult_BEST4+ epithelial",
  "Adult_Colonocyte", 
  "Adult_Goblet cell",
  "Adult_Microfold cell",
  "Adult_Tuft",
  "Fetal_Stem cells" ,
  "Fetal_TA",
  "Fetal_BEST4+ epithelial",
  "Fetal_Colonocyte","Fetal_Distal progenitor", 
  "Fetal_Enterocyte",
  "Fetal_Goblet cell", 
  "Fetal_Proximal progenitor",
  "Adult_B cells",
  "Adult_Myeloid",
  "Adult_Plasma cells",
  "Adult_T cells",
  "Fetal_B cells" , "Fetal_Red blood cells",
  "Fetal_T cells", "Fetal_Myeloid",
  "Adult_Enteroendocrine",
  "Adult_Endothelial",
  "Adult_Mesenchymal",
  "Adult_Neuronal",
  "Fetal_Endothelial",
  "Fetal_Enteroendocrine",
  "Fetal_Neuronal",
  "Fetal_Mesenchymal",
  "Cortex")))

######################
# Plot global heatmap
x = (1+exp(-pr_fullref$logitMat))**-1
#Do softmax
x = x/rowSums(x)
#Collapse column groupings into one column unless specified below
s = as.character(gut_colsplit_df$wide_annot[match(colnames(x),gut_colsplit_df$celltype)])
#Create a second level of categories that splits the adult epithelium as per the cell signal analysis
#Stem and TA
w = which(as.character(gut_colsplit_df$celltype) %in% c("Adult_Stem.cells", "Adult_TA"))
s[w] = 'Adult_stemEpithelium'
#Specialised and mature epithelium
w = which(as.character(gut_colsplit_df$celltype) %in% c('Adult_BEST4..epithelial','Adult_Colonocyte','Adult_Tuft','Adult_Goblet.cell','Adult_Paneth','Adult_Microfold.cell','Adult_BEST2..Goblet.cell'))
s[w] = 'Adult_otherEpithelium'
s = split(colnames(x),s)
x = ct$collapseClasses(x,setNames(rep(names(s),lengths(s)),unlist(s))[colnames(x)],collapseFun=sum)
#Order in the way we'd like it shown
x = x[,c('Early embryo','Gastrulation',#The precursor block
         'Fetal_epithelial','Adult_stemEpithelium','Adult_otherEpithelium',#The epithelial block
         'Fetal_Fibroblasts','Adult_Fibroblasts',#Fibro block
         'Fetal_Endothelial','Adult_Endothelial',#Endothelial block
         'Fetal_leukocytes','Adult_leukocytes','Other blood cells',#Immune block
         'Cortex')]
#Define how to group within this order
colSplits = rep(c('EarlyEmbryo','Epithelial','Fibroblasts','Endothelium','Hematopoeitic','Neg. Control'),c(2,3,2,2,3,1))
colSplits = factor(colSplits,levels=unique(colSplits))
#Do the plotting
hm = ct$plotByGroup(x,factor(crc_lee_srat@meta.data[rownames(x),'annotPlot'],levels=c(unique(grep('_SMC[0-9]+',crc_lee_srat@meta.data$annotPlot,value=TRUE)),'Endothelium','Leukocytes','Stromal cells')),
                 cluster_rows=FALSE,
                 cluster_columns=FALSE,
                 column_title_rot=90,
                 column_split=colSplits,
                 )
pdf(file.path(plotDir,"gut_lr.pdf"), height = 5, width=6,useDingbats = F)
draw(hm)
dev.off()




