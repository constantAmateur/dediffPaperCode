#' Create (some of) the figure panels relating to the analysis of the liver data.  Based on Gerda's code and data...
setwd('~/trueHome/Projects/DeDiff')

#############
# Libraries #
#############

library(SoupX)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ComplexHeatmap)
library(Matrix)
import('Code/cellSignalAnalysis.R',as='csa')
import('Code/logisticRegressionCellTypist.R',as='ct')


##########
# Params #
##########

plotDir = 'Results/LiverFigures/'
s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes
srcLiverAndCortex = 'Results/liver_and_cortex.rds'
srcLiverGastEE = 'Results/liver_gast_ee_ref.rds'
srcInHouseNB = 'Data/GOSH020'
srcHCC = 'Data/scHCC'
srcLiverCS = 'Results/cellSignalOutput/Gerdas/liverFull'
nParallel=48
fitModelPaths = c(liverAndCortex = 'Results/liverAndCortexCT.pkl',
                  liverGastEE = 'Results/liverGastEE_CT.pkl',
                  liverGastEE_HCC = 'Results/liverGastEE_HCC_CT.pkl')


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

########################
# Cell signal analysis #
########################

plotType='boxplotsByExposure'
#Load exposures.  Gerda ran these separately because of her obsession with not loosing too many transcripts in the intercept operation.  So the second is just the CSA of the 2 fetal samples.
liver_cellsig=csa$normaliseExposures(file.path(srcLiverCS,'allLiver_OutRun_fitExposures.tsv'))
fetal_liver_cellsig=csa$normaliseExposures(file.path(srcLiverCS,'fetalLiver_OutRun_fitExposures.tsv'))
#Get metadata
liver_meta = read.table(file.path(srcLiverCS,'liver_bulk_sample_metadata.txt'),sep='\t',header=TRUE)
rownames(liver_meta)=liver_meta$sample_name_no_dashes
liver_meta_rel=liver_meta[colnames(liver_cellsig$exposures),]
#Make supplementary plot
pdf(file.path(plotDir,'cellSignalLiverFull.pdf'),width=9,height=7)
hm = csa$plotExposures(liver_cellsig, 
                       column_split=liver_meta_rel$dataset_and_sample_type,
                       show_column_names=F, 
                       show_column_dend=F,
                       column_title_rot=90, 
                       cluster_column_slices=F,
                       height = nrow(liver_cellsig$exposures)*unit(2, "mm"), column_gap=unit(3, "mm"))
draw(hm)
dev.off()
liver_cellsig_df=as.data.frame(liver_cellsig$exposures)
fetal_liver_cellsig_df=as.data.frame(fetal_liver_cellsig$exposures)
#Combine them
lc_df=as.data.frame(t(cbind(liver_cellsig_df, fetal_liver_cellsig_df[rownames(liver_cellsig_df),])))
#Make summary columns
lc_df$total_adult_signal=rowSums(lc_df[,grep("Adult", colnames(lc_df))])
lc_df$total_fetal_signal=rowSums(lc_df[,grep("Fetal", colnames(lc_df))])
lc_df$earlyEmbryo = rowSums(lc_df[,c("Epiblast","Hypoblast")])
lc_df$gastrulation = rowSums(lc_df[,c("Endoderm","Mesoderm","Non-Neural Ectoderm")])
lc_df$otherEpithelium = rowSums(lc_df[,c('Fetal_HHyP','Adult_HHyP','Adult_BECs')])
#Add in the meta-data
lc_df$dataset=liver_meta$dataset
lc_df$sample_type=liver_meta$sample_type
lc_df=lc_df[lc_df$sample_type%in%c("Hepatocellular Carcinoma", "Solid Tissue Normal",
                                               "Hepatoblastoma","GTEx (brain)","Hanley (fetal liver)"),]
lc_df$sample=rownames(lc_df)
#Final cleanup of namings
lc_df$sample_type[which(lc_df$sample_type=="Hanley (fetal liver)")] = "Normal liver (fetal)"
lc_df$sample_type[which(lc_df$sample_type=="Solid Tissue Normal")] = "Normal liver (adult)"
lc_df$sample_type=factor(lc_df$sample_type, levels= c("GTEx (brain)", "Normal liver (fetal)",  "Normal liver (adult)",
                                                                            "Hepatoblastoma", "Hepatocellular Carcinoma"))
noms = colnames(lc_df)
noms[noms=="Intercept"] = "Unexplained signal"
noms[noms=="Cortex"] = "Negative control\n(adrenal cortex)"
noms[noms=="total_fetal_signal"] = "Total\nfetal signal"
noms[noms=="total_adult_signal"] = "Total\nadult signal"
noms[noms=="Fetal_Hepatocyte"] = "Fetal\nhepatocyte"
noms[noms=="Adult_Hepatocytes"] = "Adult\nhepatocyte"
noms[noms=="Unexplained signal"] = "Unexplained\nsignal"
noms[noms=="earlyEmbryo"] = "Early Embryo\nsignal"
noms[noms=="gastrulation"] = "Gastrulation\nsignal"
noms[noms=='otherEpithelium'] = 'Other Epithelial\nsignal'
colnames(lc_df) = noms
#Final ordering and pull-out of the ones we care about
dd = lc_df[,c("Early Embryo\nsignal","Gastrulation\nsignal","Fetal\nhepatocyte","Adult\nhepatocyte",'Other Epithelial\nsignal',"Unexplained\nsignal","Negative control\n(adrenal cortex)",'sample_type')]
############
# Make plot
tmp = table(dd$sample_type)
#To generate the adult colours, use colorRampPalette around this colour '#CC79A7'
colfunc <- colorRampPalette(c('#E0C6D5', "#7F1153"))
colMap = c("Unexplained\nsignal"='lightgrey',
           "Early Embryo\nsignal"='#4DBBD5',
           "Gastrulation\nsignal"='#8F3931',
           "Fetal\nhepatocyte"='#8888FF',
           "Adult\nhepatocyte"='#CC79A7',
           'Other Epithelial\nsignal'="#7F1153",
           "Negative control\n(adrenal cortex)"='black')
pdf(file.path(plotDir,'cellSignalLiver.pdf'),width=7,height=3.5)
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

########################
# Correlation with AFP #
########################

mDat_liver=readRDS(file.path(srcLiverCS,'mDat_liver.rds'))
liver_cellsig_df=as.data.frame(t(liver_cellsig$exposures))
liver_cellsig_df$best_lab=colnames(liver_cellsig_df)[apply(liver_cellsig_df,1,which.max)]
liver_cellsig_tcga=liver_cellsig_df[1:424,]
liver_cellsig_tcga$sample_name=rownames(liver_cellsig_tcga)
liver_cellsig_tcga$tcga_meta_sample_name=rownames(mDat_liver)
liver_cellsig_tcga$sample_type=mDat_liver$cgc_sample_sample_type
#get the bulk liver counts and convert to TPMs, then add to the liver_cellsig_TCGA file
rc2 = readRDS('Data/rse_gene_TCGA.RDS')
liver_se=rc2[,rownames(mDat_liver)]
gene_info=data.frame( length=liver_se@rowRanges$bp_length, kb_length=liver_se@rowRanges$bp_length/1000)
count_mtx=liver_se@assays$data$counts
count_mtx_rpk=apply(count_mtx, 2, function(x){x/gene_info$kb_length})
scale_factor=colSums(count_mtx_rpk)/1000000
count_mtx_tpm=t(apply(count_mtx_rpk, 1, function(x){x/scale_factor}))
#grep("ENSG00000167701",rownames(count_mtx_tpm))
afp_and_liver_function_genes=c("ENSG00000081051.7","ENSG00000167701.13", "ENSG00000120053.10", "ENSG00000162551.13", "ENSG00000100031.18")
liver_cellsig_tcga$AFP_tpm=count_mtx_tpm["ENSG00000081051.7",]
liver_cellsig_tcga$ALT_tpm=count_mtx_tpm["ENSG00000167701.13",]
liver_cellsig_tcga$AST_tpm=count_mtx_tpm["ENSG00000120053.10",]
liver_cellsig_tcga$ALP_tpm=count_mtx_tpm["ENSG00000162551.13",]
liver_cellsig_tcga$GGT_tpm=count_mtx_tpm["ENSG00000100031.18",]
just_tums=mDat_liver
rownames(just_tums)=rownames(liver_cellsig_tcga)
liver_cellsig_tum=liver_cellsig_tcga[which(liver_cellsig_tcga$sample_type%in%c( "Primary Tumor","Recurrent Tumor")),]
tcga_liver_tum_meta=just_tums[rownames(liver_cellsig_tum),]
meta_df=as.data.frame(tcga_liver_tum_meta)
rel_meta=data.frame(sample=rownames(liver_cellsig_tum),
                    stage=tcga_liver_tum_meta$xml_stage_event_pathologic_stage, 
                    location=tcga_liver_tum_meta$xml_tumor_tissue_site,
                    vital_status=tcga_liver_tum_meta$xml_vital_status, 
                    disease_type=tcga_liver_tum_meta$cgc_file_disease_type,
                    percent_tum_nuc=tcga_liver_tum_meta$cgc_slide_percent_tumor_nuclei,
                    percent_tum_cell=tcga_liver_tum_meta$cgc_slide_percent_tumor_cells,
                    percent_stromal=tcga_liver_tum_meta$cgc_slide_percent_stromal_cells,
                    percent_normal=tcga_liver_tum_meta$cgc_slide_percent_normal_cells,
                    percent_necrosis=tcga_liver_tum_meta$cgc_slide_percent_necrosis,
                    pathologic_t=tcga_liver_tum_meta$cgc_case_pathologic_t,
                    pathologic_n=tcga_liver_tum_meta$cgc_case_pathologic_n,
                    gender=tcga_liver_tum_meta$xml_gender,
                    age_at_diagnosis=tcga_liver_tum_meta$xml_age_at_initial_pathologic_diagnosis,
                    days_to_death=tcga_liver_tum_meta$xml_days_to_death,
                    case_tumor_status=tcga_liver_tum_meta$cgc_case_tumor_status,
                    follow_up_tumour_status=tcga_liver_tum_meta$cgc_follow_up_tumor_status,
                    follow_up_vital_status=tcga_liver_tum_meta$cgc_follow_up_vital_status,
                    case_histological_diagnosis=tcga_liver_tum_meta$cgc_case_histological_diagnosis,
                    new_tum_after_treatment=tcga_liver_tum_meta$cgc_follow_up_new_tumor_event_after_initial_treatment,
                    neopl_status=tcga_liver_tum_meta$xml_person_neoplasm_cancer_status, 
                    has_new_tum_events=tcga_liver_tum_meta$xml_has_new_tumor_events_information,
                    ethnicity=tcga_liver_tum_meta$xml_race_list, 
                    bmi=tcga_liver_tum_meta$gdc_cases.exposures.bmi,
                    histologic_grade=tcga_liver_tum_meta$xml_neoplasm_histologic_grade,
                    family_cancer_history=tcga_liver_tum_meta$xml_relative_family_cancer_history, 
                    risk_factors=tcga_liver_tum_meta$xml_history_hepato_carcinoma_risk_factors, 
                    fetoprotein_outcome_value=tcga_liver_tum_meta$xml_fetoprotein_outcome_value,
                    fibrosis_score=tcga_liver_tum_meta$xml_fibrosis_ishak_score, 
                    adjacent_inflamation=tcga_liver_tum_meta$xml_adjacent_hepatic_tissue_inflammation_extent_type,
                    viral_serologies=tcga_liver_tum_meta$xml_viral_hepatitis_serologies, 
                    fetal_hep=liver_cellsig_tum$Fetal_Hepatocyte,
                    adult_hep=liver_cellsig_tum$Adult_Hepatocytes,
                    adult_hhyp=liver_cellsig_tum$Adult_HHyP,
                    fetal_hhyp=liver_cellsig_tum$Fetal_HHyP,
                    intercept=liver_cellsig_tum$Intercept, 
                    gast_sig=liver_cellsig_tum$Endoderm+liver_cellsig_tum$Mesoderm+liver_cellsig_tum$`Non-Neural Ectoderm`,
                    best_lab=liver_cellsig_tum$best_lab, 
                    AFP_tpm=liver_cellsig_tum$AFP_tpm,
                    ALT_tpm=liver_cellsig_tum$ALT_tpm,
                    AST_tpm=liver_cellsig_tum$AST_tpm,
                    ALP_tpm=liver_cellsig_tum$ALP_tpm,
                    GGT_tpm=liver_cellsig_tum$GGT_tpm)
#Gerda's one is good enough, nothing to tweak here
pdf(file.path(plotDir,'AFP.pdf'),width=4,height=3.5)
gg = ggplot(rel_meta, aes(y=fetal_hep, x=log(AFP_tpm))) +geom_point(shape=19, alpha=0.2) +
  geom_smooth(method='lm', formula= y~x, color="#001c55") +
  stat_cor(label.y = 0.6, aes(label =  ..rr.label..))+ #this means at 35th unit in the y axis, the r squared and p value will be shown
  stat_regline_equation(label.y = 0.5) + #this means at 30th unit regresion line equation will be shown
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="AFP expression (log(TPM))", y="Fetal hepatocyte signal")
plot(gg)
dev.off()
#summary(lm(log(AFP_tpm+1)~fetal_hep, data=rel_meta))



##################################
# Load pre-calculated references #
##################################

#These are all created as per Gerda's 6_liver_single_cell.r code.
#liver_and_cortex = readRDS(srcLiverAndCortex)
liver_gastr_ee = readRDS(srcLiverGastEE)


##########################
# InHouse HBB processing #
##########################

#Get the inhouse HB data
paths = file.path(list.files(srcInHouseNB,full.names=TRUE),'filtered_feature_bc_matrix')
inhouse_hb_mtx=Read10X(paths)
inhouse_hb_srat=CreateSeuratObject(inhouse_hb_mtx)
inhouse_hb_srat[["percent.mt"]] = PercentageFeatureSet(inhouse_hb_srat, pattern = "^MT-")
#VlnPlot(inhouse_hb_srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
process_seurat =function(srat){
  srat[["percent.mt"]] = PercentageFeatureSet(srat, pattern = "^MT-")
  srat=subset(srat, subset = nFeature_RNA > 200 & nCount_RNA > 500  & percent.mt < 20)
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
inhouse_hb_srat=process_seurat(inhouse_hb_srat)
lit_genes=c('PTPRC','PDGFRB','EPCAM','PECAM1','HBA1','HBB')
#annotate the cells based on marker expression and automated annot
inhouse_hb_srat@meta.data$annot="x"
inhouse_hb_srat@meta.data$annot[which(inhouse_hb_srat@meta.data$seurat_clusters%in%c(8,9,10,11,13,14,18,19,20,21,22))]="Leukocytes"
inhouse_hb_srat@meta.data$annot[which(inhouse_hb_srat@meta.data$seurat_clusters%in%c(12))]="Erythroid cells"
inhouse_hb_srat@meta.data$annot[which(inhouse_hb_srat@meta.data$seurat_clusters%in%c(4))]="LSECs"
inhouse_hb_srat@meta.data$annot[which(inhouse_hb_srat@meta.data$seurat_clusters%in%c(15))]="Hepatic stellate cells"
inhouse_hb_srat@meta.data$annot[which(inhouse_hb_srat@meta.data$seurat_clusters%in%c(0,1,2,3,5,6,7,16,17))]="Cancer cells"
inhouse_hb_srat@meta.data$annot[which(inhouse_hb_srat@meta.data$seurat_clusters%in%c(23))]="Unknown"
#DimPlot(inhouse_hb_srat, group.by = "annot")
inhouse_hb_srat@meta.data$annot=factor(inhouse_hb_srat@meta.data$annot,
                                     levels=c("Cancer cells","Hepatic stellate cells", 
                                              "LSECs", "Erythroid cells",
                                              "Leukocytes","Unknown"))
colfunc <- colorRampPalette(c("#084594", "#C6DBEF")) 
inhouse_hb_srat@meta.data$plotAnnot = as.character(inhouse_hb_srat@meta.data$annot)
inhouse_hb_srat@meta.data$plotAnnot[which(inhouse_hb_srat@meta.data$plotAnnot %in% c('Hepatic stellate cells'))] = 'Stromal cells'
inhouse_hb_srat@meta.data$plotAnnot[which(inhouse_hb_srat@meta.data$plotAnnot %in% c('LSECs'))] = 'Endothelium'
#Plot time!
pdf(file.path(plotDir,"hb_inhouse_umap.pdf"), height = 3, width=3,useDingbats = F)
gg=DimPlot(inhouse_hb_srat, group.by = "plotAnnot", label=T,
         #cols = c("#084594","#7E481C","#CC6677" , "#771122","#DDCC77", "lightgrey"),
         cols = c("#084594",'#A5636D',"#501D28","#DDCC77","#44AA99", "lightgrey"),
         pt.size = 0.1) +theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
                               axis.ticks = element_blank(), legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1") +ggtitle("")
plot(gg)
dev.off()


##################
# HCC Processing #
##################

hcc_sc_paths = file.path(grep('SRR[0-9]+$',list.files(srcHCC,full.names=TRUE),value=TRUE),'output','Gene','filtered')
hcc_sc_mtx=Read10X(hcc_sc_paths)
hcc_sc_srat=CreateSeuratObject(hcc_sc_mtx)
hcc_sc_srat@meta.data$sample_no=Idents(hcc_sc_srat)
process_seurat =function(srat){
  srat[["percent.mt"]] = PercentageFeatureSet(srat, pattern = "^MT-")
  srat=subset(srat, subset = nFeature_RNA > 300 & nCount_RNA > 1000  & percent.mt < 10)
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
hcc_sc_srat=process_seurat(hcc_sc_srat)
#DimPlot(hcc_sc_srat, label=T)
#annotate based on markers from the Ho et al. paper (since annotation not provided)
hcc_sc_srat@meta.data$annot="cancer"
hcc_sc_srat@meta.data$annot[which(hcc_sc_srat@meta.data$seurat_clusters%in%c(3,4,5,8,9,10,11,12,13,14,
                                                                             16,18,20,21,24,25, 27,30,34))]="Leukocytes"
hcc_sc_srat@meta.data$annot[which(hcc_sc_srat@meta.data$seurat_clusters%in%c(36))]="Endothelium"
hcc_sc_srat@meta.data$annot2=hcc_sc_srat@meta.data$annot
hcc_sc_srat@meta.data$annot2[which(hcc_sc_srat@meta.data$annot=="cancer")]=paste0(hcc_sc_srat@meta.data$annot[which(hcc_sc_srat@meta.data$annot=="cancer")], "_",
                                                                                  hcc_sc_srat@meta.data$seurat_clusters[which(hcc_sc_srat@meta.data$annot=="cancer")])
hcc_sc_srat@meta.data$annot3=hcc_sc_srat@meta.data$annot
hcc_sc_srat@meta.data$annot3[which(hcc_sc_srat@meta.data$annot=="cancer")]=paste0(hcc_sc_srat@meta.data$annot[which(hcc_sc_srat@meta.data$annot=="cancer")], "_",
                                                                                  hcc_sc_srat@meta.data$orig.ident[which(hcc_sc_srat@meta.data$annot=="cancer")])
colfunc <- colorRampPalette(c("#084594", "#C6DBEF"))
#Plot time!
pdf(file.path(plotDir,"hcc_sc_umap.pdf"), height = 3, width=3,useDingbats = F)
gg=DimPlot(hcc_sc_srat, group.by = "annot3", label=T,
           #cols = c( colfunc(8),"#771122","#DDCC77"),
           cols = c( colfunc(8),'#A5636D',"#DDCC77"),
           pt.size = 0.1) +theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
                                 axis.ticks = element_blank(), legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1") +ggtitle("")
plot(gg)
dev.off()

#######################
# Logistic Regression #
#######################

#############
# Fit models
#Define shared gene-sets for LR
inhouse_ee_gast_shared_genes=intersect(rownames(inhouse_hb_srat), rownames(liver_gastr_ee))
hcc_inter_ee_shared_genes = intersect(rownames(hcc_sc_srat), rownames(liver_gastr_ee))
#Train the models and fit
#Inhouse
tgtFile = fitModelPaths['liverGastEE']
if(!file.exists(tgtFile)){
  ct$trainCelltypistModel(liver_gastr_ee@assays$RNA@counts[inhouse_ee_gast_shared_genes, ],
                          outPath = tgtFile,
                          n_jobs = nParallel,
                          labels=as.character(liver_gastr_ee@meta.data$annot))
}
pr_inh_ee = ct$runCelltypist(inhouse_hb_srat@assays$RNA@counts[inhouse_ee_gast_shared_genes,],
                          model=tgtFile)
#HCC
tgtFile = fitModelPaths['liverGastEE_HCC']
if(!file.exists(tgtFile)){
  ct$trainCelltypistModel(liver_gastr_ee@assays$RNA@counts[hcc_inter_ee_shared_genes, ],
                          outPath = tgtFile,
                          n_jobs = nParallel,
                          labels=as.character(liver_gastr_ee@meta.data$annot))
}
pr_hcc_ee = ct$runCelltypist(hcc_sc_srat@assays$RNA@counts[hcc_inter_ee_shared_genes,],
                             model=fitModelPaths['liverGastEE_HCC'])

##########################
# Define column groupings
x = pr_inh_ee$logitMat
colsplit_df=data.frame(celltype=colnames(x))
colsplit_df$wide_annot=colsplit_df$celltype
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Fetal_Mono-Mac","Fetal_B-cells","Fetal_DCs",
                                                         "Fetal_Monocyte",
                                                         "Fetal_Early lymphoid_T lymphocyte",
                                                         "Fetal_VCAM1+ EI macrophage", "Fetal_ILC precursor",
                                                         "Fetal_Neutrophil-myeloid progenitor",
                                                         "Fetal_Monocyte precursor", 
                                                         "Fetal_NK",
                                                         "Fetal_Mast cell", 
                                                         "Fetal_Kupffer Cell")))]="Fetal_leukocytes"
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Adult_Non-inflammatory_macs", 
                                                         "Adult_T-cells", 
                                                         "Adult_Inflammatory_macs",
                                                         "Adult_NK-like_cells",
                                                         "Adult_Mature_Bcells")))]="Adult_leukocytes"
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Adult_Hepatocytes", 
                                                         "Fetal_Hepatocyte", 
                                                         "Adult_HHyP",
                                                         "Fetal_HHyP",
                                                         "Adult_BECs")))]="EPCAM+"
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Fetal_Erythroid", 
                                                         "Fetal_MEMP", 
                                                         "Fetal_Megakaryocyte",
                                                         "Fetal_MEMP",
                                                         "Fetal_HSC_MPP",
                                                         "Adult_Plasma_cells",
                                                         "Adult_Erythroid_cells")))]="Other blood cells"
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Cortex")))]="Negative control"
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Hypoblast", "Epiblast")))]="Early embryo"
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Mesoderm", "Endoderm", "Non-Neural Ectoderm")))]="Gastrulation"
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Fetal_Endothelial cell")))]="Fetal_Endothelium"
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%gsub('[ -]|\\+','.',c("Adult_LSECs",
                                                         "Adult_Portal_endothelium")))]="Adult_Endothelium"
colsplit_df$wide_annot=factor(colsplit_df$wide_annot, levels=c("Early embryo", "Gastrulation","EPCAM+", "Adult_Hepatic_stellate_cells", "Fetal_Endothelium",'Adult_Endothelium',
                                                               "Other blood cells","Adult_leukocytes", "Fetal_leukocytes", 
                                                               "Fetal_Fibroblast",  "Negative control" ))
colsplit_df$celltype=factor(colsplit_df$celltype, levels =gsub('[ -]|\\+','.', c("Hypoblast", "Epiblast","Mesoderm", "Endoderm", "Non-Neural Ectoderm","Fetal_Hepatocyte",
                                                             "Adult_Hepatocytes",
                                                             "Fetal_HHyP",
                                                             "Adult_HHyP",
                                                             "Adult_BECs", 
                                                             "Adult_Hepatic_stellate_cells", 
                                                             "Fetal_Endothelial cell",
                                                             "Adult_LSECs",
                                                             "Adult_Portal_endothelium",
                                                             "Adult_Erythroid_cells",
                                                             "Fetal_Erythroid", 
                                                             "Fetal_MEMP", 
                                                             "Fetal_Megakaryocyte",
                                                             "Fetal_HSC_MPP",
                                                             "Adult_Plasma_cells",
                                                             "Adult_Non-inflammatory_macs", 
                                                             "Adult_T-cells", 
                                                             "Adult_Inflammatory_macs",
                                                             "Adult_NK-like_cells",
                                                             "Adult_Mature_Bcells",
                                                             "Fetal_Mono-Mac","Fetal_B-cells","Fetal_DCs",
                                                             "Fetal_Monocyte",
                                                             "Fetal_Early lymphoid_T lymphocyte",
                                                             "Fetal_VCAM1+ EI macrophage", "Fetal_ILC precursor",
                                                             "Fetal_Neutrophil-myeloid progenitor",
                                                             "Fetal_Monocyte precursor", 
                                                             "Fetal_NK",
                                                             "Fetal_Mast cell", 
                                                             "Fetal_Kupffer Cell",
                                                             "Fetal_Fibroblast",
                                                             "Cortex")))

######################
# Plot global heatmap
#Get the things worth plotting
x1 = pr_inh_ee$logitMat
#x1 = (1+exp(-pr_inh_ee$logitMat))**-1
x1 = x1[!inhouse_hb_srat@meta.data[rownames(x1),'annot'] %in% c('Unknown'),]
x2 = pr_hcc_ee$logitMat
#x2 = (1+exp(-pr_hcc_ee$logitMat))**-1
#Make immune label disease specific
noms1 = as.character(inhouse_hb_srat@meta.data[rownames(x1),'annot'])
noms1[noms1=='Leukocytes'] = 'HB_Leukocytes'
noms2 = as.character(hcc_sc_srat@meta.data[rownames(x2),'annot3'])
noms2[noms2=='Leukocytes'] = 'HCC_Leukocytes'
mDat = data.frame(cellName=paste0(rep(c('HB_','HCC_'),c(nrow(x1),nrow(x2))),c(rownames(x1),rownames(x2))),
                  cellType = c(noms1,noms2),
                  disease = rep(c('HB','HCC'),c(nrow(x1),nrow(x2))))
x = rbind(x1,x2)
#Prob space plus softmax
x = (1+exp(-x))**-1
x = x/rowSums(x)
#Simplify columns
s = as.character(colsplit_df$wide_annot[match(colnames(x),colsplit_df$celltype)])
s = ifelse(s %in% c("Adult_Hepatic_stellate_cells", "Fetal_Fibroblast", "Negative control", "EPCAM+"),colnames(x),s)
#Add in grouping of HHyP and BECs into one category
s[s %in% c('Adult_HHyP','Fetal_HHyP','Adult_BECs')] = 'Other_Epithelial'
s = split(colnames(x),s)
x = ct$collapseClasses(x,setNames(rep(names(s),lengths(s)),unlist(s))[colnames(x)],collapseFun=sum)
#x = do.call(cbind,lapply(s,function(e) apply(x[,e,drop=FALSE],1,max)))
#Order in the way we'd like it shown
x = x[,c('Early embryo','Gastrulation',#The precursor block
         # 'Fetal_Hepatocyte','Adult_Hepatocytes','Fetal_HHyP','Adult_HHyP','Adult_BECs',#The epithelial block
         'Fetal_Hepatocyte','Adult_Hepatocytes','Other_Epithelial',#The epithelial block
         'Fetal_Fibroblast','Adult_Hepatic_stellate_cells',#Fibro block
         'Fetal_Endothelium','Adult_Endothelium',#Endothelial block
         'Fetal_leukocytes','Adult_leukocytes','Other blood cells', #Immune block
         'Cortex')]
#Define how to group within this order
#colSplits = rep(c('EarlyEmbryo','Epithelial','Fibroblasts','Endothelium','Hematopoeitic','Neg. Control'),c(2,5,2,2,3,1))
colSplits = rep(c('EarlyEmbryo','Epithelial','Fibroblasts','Endothelium','Hematopoeitic','Neg. Control'),c(2,3,2,2,3,1))
colSplits = factor(colSplits,levels=unique(colSplits))
#Do the plotting
hm = ct$plotByGroup(x,mDat$cellType,
                 cluster_rows=FALSE,
                 cluster_columns=FALSE,
                 column_title_rot=90,
                 column_split=colSplits,
                 row_split = ifelse(sort(unique(mDat$cellType)) %in% noms1,'HB','HCC'))
pdf(file.path(plotDir,"liver_lr.pdf"), height = 5, width=6,useDingbats = F)
draw(hm)
dev.off()

