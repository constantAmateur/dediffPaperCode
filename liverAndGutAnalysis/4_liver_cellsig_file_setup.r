source('process_seurat.R')
library(SummarizedExperiment)
library(ggplot2)
library(reshape2)
library(stringr)

##############Fetal liver reference - from Muzz's lab##########################
l_f_mtx=readMM('/lustre/scratch117/casm/team274/gk14/Dediff/fetal_liver/matrix.mtx')
l_f_genes=read.csv('/lustre/scratch117/casm/team274/gk14/Dediff/fetal_liver/obs.csv')
l_f_meta=read.csv('/lustre/scratch117/casm/team274/gk14/Dediff/fetal_liver/var.csv')
rownames(l_f_meta)=l_f_meta$index

rownames(l_f_mtx)=l_f_genes$index
colnames(l_f_mtx)=l_f_meta$index

fetal_liver_srat=CreateSeuratObject(l_f_mtx, meta.data = l_f_meta)

fetal_liver_srat_processed=process_seurat(fetal_liver_srat)

DimPlot(fetal_liver_srat_processed, group.by = "cell.labels")

######################adult liver - GSE115469###################################
adult_liver_paths=c('/lustre/scratch117/casm/team274/gk14/Dediff/GSE115469_adult_liver/P1TLH/',
                    '/lustre/scratch117/casm/team274/gk14/Dediff/GSE115469_adult_liver/P2TLH/',
                    '/lustre/scratch117/casm/team274/gk14/Dediff/GSE115469_adult_liver/P3TLH/',
                    '/lustre/scratch117/casm/team274/gk14/Dediff/GSE115469_adult_liver/P4TLH/',
                    '/lustre/scratch117/casm/team274/gk14/Dediff/GSE115469_adult_liver/P5TLH/')
names(adult_liver_paths)=c("P1TLH", "P2TLH","P3TLH","P4TLH","P5TLH")

adult_liver=Read10X(adult_liver_paths)
colnames(adult_liver)=str_replace(string = colnames(adult_liver), pattern = "-", replacement = "_")
adult_liver_clusts=read.table('/lustre/scratch117/casm/team274/gk14/Dediff/GSE115469_adult_liver/Cell_clusterID_cycle.txt', sep = "\t", header = T)
rownames(adult_liver_clusts)=adult_liver_clusts$CellName
adult_liver=adult_liver[,rownames(adult_liver_clusts)]
adult_liver_srat=CreateSeuratObject(adult_liver, meta.data = adult_liver_clusts)

adult_liver_srat_processed=process_seurat(adult_liver_srat)

DimPlot(adult_liver_srat_processed, group.by = "annot")

clusts = 1:20
annot= c("Hep_1","ab_Tcells","Hep_2","Inflammatory_macs", "Hep_3",
         "Hep_4", "Plasma_cells", "NK-like_cells", "gd_Tcells_1", "Non-inflammatory_macs", 
         "Periportal_LSECs", "Central_venous_LSECs", "Portal_endothelium", "Hep_5", "Hep_6",
         "Mature_Bcells", "Cholangiocytes", "gd_Tcells_2", "Erythroid_cells", "Hepatic_stellate_cells")

adult_liver_srat@meta.data$Cluster.[adult_liver_srat@meta.data$Cluster. %in% clusts] = annot[match(adult_liver_srat@meta.data$Cluster., clusts, nomatch = 0)]

adult_liver_srat_processed$annot=adult_liver_srat@meta.data$Cluster.



#merge some of adult reference cell labels 

adult_liver_srat_processed@meta.data$annot[grep("Hep_", adult_liver_srat_processed@meta.data$annot)]="Hepatocytes"
adult_liver_srat_processed@meta.data$annot[grep("Tcells", adult_liver_srat_processed@meta.data$annot)]="T-cells"
adult_liver_srat_processed@meta.data$annot[grep("LSECs", adult_liver_srat_processed@meta.data$annot)]="LSECs"
#merge some of fetal reference cell labels 
fetal_liver_srat@meta.data$annot=fetal_liver_srat@meta.data$cell.labels
fetal_liver_srat@meta.data$annot[grep("Erythroid", fetal_liver_srat@meta.data$annot)]="Erythroid"
fetal_liver_srat@meta.data$annot[grep("B cell", fetal_liver_srat@meta.data$annot)]="B-cells"
fetal_liver_srat@meta.data$annot[grep("DC", fetal_liver_srat@meta.data$annot)]="DCs"

#add fetal or adult cell type
fetal_liver_srat@meta.data$annot=paste0("Fetal_",fetal_liver_srat@meta.data$annot)
adult_liver_srat_processed@meta.data$annot=paste0("Adult_",adult_liver_srat_processed@meta.data$annot)

adr_genes=read.table("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/SCP/Data/fAdrenal19wk/cellranger302_count_32644_WSSS_F_Adr8710632_GRCh38-1_2_0/filtered_feature_bc_matrix/features.tsv.gz")

adr_fliver=adr_genes
rownames(adr_fliver)=make.unique(adr_fliver$V2)
adr_and_fliver_inter=intersect(rownames(fetal_liver_srat), adr_fliver$V2)
adr_fliver=adr_fliver[adr_and_fliver_inter,]
liver_ens_mtx_f=fetal_liver_srat@assays$RNA@counts[adr_and_fliver_inter,]
rownames(liver_ens_mtx_f)=adr_fliver$V1

adr_aliver=adr_genes
rownames(adr_aliver)=adr_aliver$V1

adult_liver_genes=read.table('/lustre/scratch117/casm/team274/gk14/Dediff/GSE115469_adult_liver/P1TLH/genes.tsv', sep = "\t")
liver_ens_mtx_a=adult_liver_srat_processed@assays$RNA@counts
rownames(liver_ens_mtx_a)=adult_liver_genes$V1
f_a_comm_ens=intersect(rownames(liver_ens_mtx_a), rownames(liver_ens_mtx_f))

new_fliver=liver_ens_mtx_f[f_a_comm_ens,]
new_aliver=liver_ens_mtx_a[f_a_comm_ens,]

f_a_comm_symbols=adr_aliver[f_a_comm_ens,]

rownames(new_fliver)=f_a_comm_symbols$V2
rownames(new_aliver)=f_a_comm_symbols$V2
new_fliver=CreateSeuratObject(new_fliver, meta.data = fetal_liver_srat@meta.data)
new_aliver=CreateSeuratObject(new_aliver, meta.data = adult_liver_srat_processed@meta.data)

merged_liver_srat=merge(new_fliver, new_aliver)


merged_liver_srat_processed=process_seurat(merged_liver_srat)
####relabel based on annotation from the hhyp paper#####

hhyp_liver=read.csv('/lustre/scratch117/casm/team274/gk14/Dediff/GSE130473_liver_hhyps/GSE130473_Series_count_matrix.csv')
rownames(hhyp_liver)=hhyp_liver$X
hhyp_liver=hhyp_liver[,2:ncol(hhyp_liver)]

hhyp_liver=as.matrix(hhyp_liver)
hhyp_liver=Matrix(hhyp_liver, sparse=T)

#get a bulk file example with gene lengths
bulk_1=read.table('matts-bits-of-code-master/cellSignalAnalysisV2/bulkData/SangerProjectTARGET_ALL_phase1/SangerProjectTARGET_ALL_phase1_SRR452447.tsv', header = T, sep = "\t")
rownames(bulk_1)=bulk_1$geneName
rownames(hhyp_liver) = sapply(strsplit(rownames(hhyp_liver), "\\."),"[", 1)
inter=intersect(rownames(hhyp_liver), rownames(bulk_1))
ensembl_ids=read.csv('onecs/obs.csv')
rownames(ensembl_ids)=ensembl_ids$gene_ids
inter2=intersect(inter, rownames(ensembl_ids))

sum(colSums(hhyp_liver[rownames(hhyp_liver)[!rownames(hhyp_liver)%in%inter2],]))
hhyp_liver=hhyp_liver[inter2,]
gene_lengths=bulk_1[inter2,]
gene_names=ensembl_ids[inter2,]

rownames(hhyp_liver)=gene_names$index
gene_lengths$gene_lengths_kb=gene_lengths$geneLengths/1000
hhyp_liver_rpk=apply(hhyp_liver, 2, function(x){x/gene_lengths$gene_lengths_kb})
hhyp_liver_rpk=Matrix(hhyp_liver_rpk, sparse = T)
hhyp_seurat=CreateSeuratObject(hhyp_liver_rpk)

process_seurat =function(srat){
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat)
  srat = RunPCA(srat, npcs = 15)
  srat = FindNeighbors(srat, dims=1:15)
  srat = FindClusters(srat, resolution = 0.6)
  srat = RunUMAP(srat, dims=1:15, min.dist = 0.5, n.neighbors = 10)
  return(srat)
}

hhyp_seurat=process_seurat(hhyp_seurat)
cell_info=read.csv('/lustre/scratch117/casm/team274/gk14/Dediff/GSE130473_liver_hhyps/GSE130473_Series_feature_matrix.csv')
cell_info$new_cell_id=cell_info$CELL_ID
cell_info$new_cell_id=str_replace(cell_info$new_cell_id, "#", ".")
rownames(cell_info)=cell_info$new_cell_id
cell_info=cell_info[colnames(hhyp_seurat),]

hhyp_seurat@meta.data$fetal_or_adult=cell_info$TISSUE
DimPlot(hhyp_seurat, group.by = "fetal_or_adult")
DimPlot(hhyp_seurat)
FeaturePlot(hhyp_seurat, c('SPP1', 'KRT7', 'TACSTD2',"CLDN4"))
DimPlot(hhyp_seurat)
FeaturePlot(hhyp_seurat, c('ALB', "AHSG", "APOA2"))

hhyp_seurat@meta.data$new_labs=as.character(hhyp_seurat@meta.data$RNA_snn_res.0.6)
#annotate the hhyp data based on markers in the paper
hhyp_seurat@meta.data[names(which(hhyp_seurat@assays$RNA@counts["MUC5B",]>5)),]$new_labs="adult_BECs"
hhyp_seurat@meta.data$new_labs[which(hhyp_seurat@meta.data$new_labs=="1")]="adult_hhyps"
hhyp_seurat@meta.data[rownames(hhyp_seurat@meta.data[hhyp_seurat@meta.data$new_labs=="adult_hhyps",])[which(hhyp_seurat@meta.data[hhyp_seurat@meta.data$new_labs=="adult_hhyps",]$fetal_or_adult=="FETAL")],]$new_labs="fetal_hhyps"
hhyp_seurat@meta.data$new_labs[which(hhyp_seurat@meta.data$new_labs%in%c("4", "8"))]="fetal_hepatocytes"
hhyp_seurat@meta.data$new_labs[which(hhyp_seurat@meta.data$new_labs%in%c("0"))]="adult_hepatocytes"

#logistic regression against the subset of relevant cell types in the 10X reference
source('logisticRegression.R')
workers=NULL
hep_and_chol=subset(merged_liver_srat_processed, annot%in%c("Fetal_Hepatocyte", "Adult_Hepatocytes", "Adult_Cholangiocytes"))

comm_genes=intersect(rownames(hep_and_chol), rownames(hhyp_seurat))
tr2=trainModel(hhyp_seurat@assays$RNA@counts[comm_genes,], as.character(hhyp_seurat@meta.data$new_labs))
pr_sc_tr2=predictSimilarity(tr2, hep_and_chol@assays$RNA@counts[comm_genes,])
similarityHeatmap(pr_sc_tr2, cluster_rows=T, cluster_columns=F,
                  row_split=hep_and_chol@meta.data$annot, row_title_rot=0, use_raster=F, 
                  column_order=c("fetal_hhyps", "adult_hhyps", "fetal_hepatocytes", "adult_hepatocytes",
                                 "adult_BECs", "2","3", "5", "6","7","9", "10", "11", "12"))

df1=as.data.frame(pr_sc_tr2)

df1$top_match=colnames(df1)[apply(df1,1,which.max)]
df1$old_celltype=hep_and_chol@meta.data$annot

hep_and_chol@meta.data$new_labs=df1$top_match

hep_and_chol@meta.data[rownames(df1[which(df1$old_celltype=="Fetal_Hepatocyte"&df1$top_match=="fetal_hhyps"),]),]$annot="Fetal_HHyP"
hep_and_chol@meta.data[rownames(df1[which(df1$old_celltype=="Adult_Cholangiocytes"&df1$top_match=="adult_hhyps"),]),]$annot="Adult_HHyP"
hep_and_chol@meta.data[rownames(df1[which(df1$old_celltype=="Adult_Cholangiocytes"&df1$top_match=="adult_BECs"),]),]$annot="Adult_BECs"

DimPlot(hep_and_chol, group.by = "new_labs")
merged_liver_srat@meta.data[rownames(df1[which(df1$old_celltype=="Fetal_Hepatocyte"&df1$top_match=="fetal_hhyps"),]),]$annot="Fetal_HHyP"
merged_liver_srat@meta.data[rownames(df1[which(df1$old_celltype=="Adult_Cholangiocytes"&df1$top_match=="adult_hhyps"),]),]$annot="Adult_HHyP"
merged_liver_srat@meta.data[rownames(df1[which(df1$old_celltype=="Adult_Cholangiocytes"&df1$top_match=="adult_BECs"),]),]$annot="Adult_BECs"
#subset to remove unsassigned "cholangiocytes" - 5 cells

merged_liver_srat=subset(merged_liver_srat, annot%in%c("Adult_Cholangiocytes"), invert=T)
saveRDS(merged_liver_srat, '/lustre/scratch117/casm/team274/gk14/Dediff/merged_liver/merged_liver_srat.rds')



#subset Cortex
adr_all=readRDS('/lustre/scratch117/casm/team274/gk14/cellxgene/adr_all.rds')
cortex=subset(adr_all, Annotation%in%c("Cortex"))
#get Cortex matrix
ensembl_adrenal=read.table("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/SCP/Data/fAdrenal19wk/cellranger302_count_32644_WSSS_F_Adr8710632_GRCh38-1_2_0/filtered_feature_bc_matrix/features.tsv.gz", sep = "\t", header = F)
cortex_mtx=cortex@assays$RNA@counts
rownames(cortex_mtx)=ensembl_adrenal$V1
#################################
#convert gene names to ensembl IDs
ensembl_ids=read.csv('onecs/obs.csv')
rownames(ensembl_ids)=ensembl_ids$index
inter=intersect(rownames(merged_liver_srat_processed),ensembl_ids$index)
ensembl_inter=ensembl_ids[inter,]
liver_mtx=merged_liver_srat_processed@assays$RNA@counts
liver_mtx=liver_mtx[inter,]

rownames(liver_mtx)= ensembl_inter$gene_ids

common_ensembl=intersect(rownames(liver_mtx), rownames(cortex_mtx))

liver_with_cortex=cbind(liver_mtx[common_ensembl,], cortex_mtx[common_ensembl,])
colnames(liver_with_cortex)=c(paste0(merged_liver_srat_processed@meta.data$annot,":",colnames(liver_mtx)),
                       paste0("Cortex:", colnames(cortex_mtx)))

writeMM(liver_with_cortex, "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_liver/liver_with_cortex_ref/liver.mtx")
write.table(rownames(liver_with_cortex), "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_liver/liver_with_cortex_ref/liver_rowNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(colnames(liver_with_cortex), "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_liver/liver_with_cortex_ref/liver_columnNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)

# get paths for liver cancers - TCGA

rc2 = readRDS('/lustre/scratch119/realdata/mdt1/team274/my4/bulkRNAseq/rse_gene_TCGA.RDS')
len=read.table('matts-bits-of-code-master/cellSignalAnalysisV2/bulkData/SangerProjectTARGET_ALL_phase1/SangerProjectTARGET_ALL_phase1_SRR452447.tsv', sep = "\t", header = T)
obs=read.csv("onecs/obs.csv")
mDat = colData(rc2)
mDat$tissue = as.character(mDat$xml_tumor_tissue_site)
table(mDat$tissue)
mDat_liver=mDat[which(mDat$tissue=="Liver"),]
mDat_liver$xml_tumor_type
mDat_liver=mDat_liver[,colnames(mDat_liver)[which(!colnames(mDat_liver)%in%names(which(apply(mDat_liver, 2, function(x){sum(is.na(x))})==ncol(mDat_liver))))]]
saveRDS(mDat_liver, "mDat_liver.rds")
paths=paste0("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/CellDeconvolution/Results/bulkRNAseq/TCGA/frag/TCGA_", rownames(mDat_liver), ".tsv")
write.table(paths, "cellsig_liver/liver_tcga_paths.txt", sep = "\t", row.names = F, col.names = F, quote = F)
#get a bulk file example with gene lengths
bulk_1=read.table('matts-bits-of-code-master/cellSignalAnalysisV2/bulkData/SangerProjectTARGET_ALL_phase1/SangerProjectTARGET_ALL_phase1_SRR452447.tsv', header = T, sep = "\t")
rownames(bulk_1)=bulk_1$geneName
#fetal liver bulk paths
hanley_fliver=read.table('Hanley Supplemental_Data_01_GENCODE18_Gene_ReadCounts.tsv',sep = "\t", header = T)

rownames(hanley_fliver)=make.unique(sapply(strsplit(hanley_fliver$gene_id, "\\."), "[", 1))
hanley_fliver$gene_id=make.unique(sapply(strsplit(hanley_fliver$gene_id, "\\."), "[", 1))
rownames(bulk_1)=bulk_1$geneName
inter=intersect(rownames(hanley_fliver), rownames(bulk_1))

hanley_fliver_inter=hanley_fliver[inter,]
bulk_1_inter=bulk_1[inter,]
hanley_fliver_inter$geneLengths=bulk_1_inter$geneLengths

liver_1=hanley_fliver_inter[,c("gene_id", "geneLengths", "Liver_1")]
liver_2=hanley_fliver_inter[,c("gene_id", "geneLengths", "Liver_2")]
colnames(liver_1)=c("geneName",  "geneLengths", "fetal_Liver_1" )
colnames(liver_2)=c("geneName",  "geneLengths", "fetal_Liver_2" )

write.table(liver_1, paste0('cellsig_liver/fetal_liver_files/fetal_Liver_1.tsv'), sep="\t", row.names = F, col.names = T, quote = F)
write.table(liver_2, paste0('cellsig_liver/fetal_liver_files/fetal_Liver_2.tsv'), sep="\t", row.names = F, col.names = T, quote = F)
#hepatoblastoma and HCC from st Judes

rc2 = readRDS('/lustre/scratch119/realdata/mdt1/team274/my4/bulkRNAseq/rse_gene_TCGA.RDS')
gene_info=rowData(rc2)
df1=data.frame(gene_symbol=unlist(gene_info$symbol), ens_id=names(unlist(gene_info$symbol)))
df1=na.omit(df1)
rownames(df1)=make.unique(df1$gene_symbol)


for (i in list.files('/home/jovyan/Dediff/stJudes_hepa/StJudes_hepatoblast', pattern= "SJ", full.names = T)) {
  fname=str_replace(basename(i), ".txt", "")
  tum_file=read.table(i, sep = "\t", header = F)
  rownames(tum_file)=tum_file$V1
  inter=intersect(tum_file$V1, rownames(df1))
  df1_inter=df1[inter,]
  tum_file=tum_file[inter,]
  
  tum_file$geneName=df1_inter$ens_id
  rownames(tum_file)=make.unique(tum_file$geneName)
  inter2=intersect(rownames(bulk_1), rownames(tum_file))
  
  bulk_1_inter=bulk_1[inter2,]
  tum_file_inter=tum_file[inter2,]
  
  tum_file_inter$geneLengths=bulk_1_inter$geneLengths
  
  tum_write_file=data.frame(geneName=tum_file_inter$geneName, geneLengths=tum_file_inter$geneLengths, fname=tum_file_inter$V2)
  colnames(tum_write_file)[3]=fname
  write.table(tum_write_file, paste0('/home/jovyan/Dediff/cellsig_liver/stJudes_files/', fname, ".tsv"), sep="\t", row.names = F, col.names = T, quote = F)
  
}

hepato_paths=list.files('/home/jovyan/Dediff/cellsig_liver/stJudes_files', full.names = T)
write.table(hepato_paths, '/home/jovyan/Dediff/cellsig_liver/hepatoblastoma_paths.txt', sep = "\t", row.names = F, col.names = F, quote = F)

#hepatoblastoma from Sekiguchi et al 

hb_seki=read.table('HBL50PC_NL10_count.txt', sep = "\t", header = T, row.names = 1)

hb_seki_meta=data.frame(sample=colnames(hb_seki), sample_type=colnames(hb_seki))

hb_seki_meta$sample_type[grep("NL", hb_seki_meta$sample_type)]="normal_liver"
hb_seki_meta$sample_type[grep("P", hb_seki_meta$sample_type)]="hepatoblastoma_pretreatment"
hb_seki_meta$sample_type[grep("C", hb_seki_meta$sample_type)]="hepatoblastoma_posttreatment"

ensembl_ids=read.csv('onecs/obs.csv')
rownames(ensembl_ids)=ensembl_ids$index

inter=intersect(rownames(hb_seki), ensembl_ids$index)
hb_seki_inter=hb_seki[inter,]
ensembl_ids_inter=ensembl_ids[inter,]

rownames(hb_seki_inter)=ensembl_ids_inter$gene_ids

inter2=intersect(rownames(hb_seki_inter), rownames(bulk_1))
hb_seki_inter2=hb_seki_inter[inter2,]
bulk_1_inter=bulk_1[inter2,]
hb_seki_inter2$geneLengths=bulk_1_inter$geneLengths
hb_seki_inter2$geneName=rownames(hb_seki_inter2)

for (i in hb_seki_meta$sample) {
  df=hb_seki_inter2[,c("geneName", "geneLengths", i)]
  write.table(df, paste0('/home/jovyan/Dediff/cellsig_liver/hb_seki/', i, ".tsv"), sep = "\t",row.names = F, col.names = T, quote = F )
  
  
}

paths_seki=paste0('/home/jovyan/Dediff/cellsig_liver/hb_seki/', hb_seki_meta$sample, ".tsv")
write.table(paths_seki, '/home/jovyan/Dediff/cellsig_liver/sekiguchi_hepatoblastoma_paths.txt', sep = "\t", row.names = F, col.names = F, quote = F)
