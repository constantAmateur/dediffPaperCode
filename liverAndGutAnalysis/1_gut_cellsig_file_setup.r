source('process_seurat.R')
library(SummarizedExperiment)
library(ggplot2)
library(reshape2)
#########gut ref with cortex as a negative ctrl##############
#load the mtx (converted from .h5ad)
gut_mtx=readMM('/lustre/scratch117/casm/team274/gk14/rasa_final_object_soupx/matrix.mtx')
gut_genes=read.csv('/lustre/scratch117/casm/team274/gk14/rasa_final_object_soupx/obs.csv')
gut_meta=read.csv('/lustre/scratch117/casm/team274/gk14/rasa_final_object_soupx/var.csv')
rownames(gut_meta)=gut_meta$X

rownames(gut_mtx)=gut_genes$gene_ids
colnames(gut_mtx)=gut_meta$X

#subset relevant timepoints and locations

gut_meta_age=gut_meta[which(gut_meta$Age_group%in%c("Adult", "First trim", "Second trim")),]
gut_meta_age_region=gut_meta_age[which(gut_meta_age$Region%in%c("APD", "LargeInt", "REC")),]


gut_rel_mtx=gut_mtx[,gut_meta_age_region$X]

#create Seurat object with relevant cells

gut_srat=CreateSeuratObject(gut_rel_mtx, meta.data = gut_meta_age_region)

saveRDS(gut_srat, '/lustre/scratch117/casm/team274/gk14/rasa_final_object_soupx/rasa_rel_srat.rds')

# create annotation where epithelial cells have narrow annotation, but everything else has broad annotation
gut_srat=readRDS('/lustre/scratch117/casm/team274/gk14/rasa_final_object_soupx/rasa_rel_srat.rds')
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
table(gut_srat$final_annot)
gut_srat=subset(gut_srat, final_annot%in%c("Fetal_CLDN10+ cells","Fetal_Microfold cell", "Fetal_Paneth", "Fetal_Tuft"), invert=T)
#subset cortex from adrenal gland dataset
adr_all=readRDS('/lustre/scratch117/casm/team274/gk14/cellxgene/adr_all.rds')
cortex=subset(adr_all, Annotation%in%c("Cortex"))
#get cortex matrix, change gene names to ensembl IDs
ensembl_adrenal=read.table("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/SCP/Data/fAdrenal19wk/cellranger302_count_32644_WSSS_F_Adr8710632_GRCh38-1_2_0/filtered_feature_bc_matrix/features.tsv.gz", sep = "\t", header = F)
cortex_mtx=cortex@assays$RNA@counts
rownames(cortex_mtx)=ensembl_adrenal$V1
#get gut matrix, change gene names to ensembl IDs
ensembl_gut=read.csv('onecs/obs.csv')
gut_rel_mtx=gut_srat@assays$RNA@counts
rownames(gut_rel_mtx)=ensembl_gut$gene_ids

comm_genes=intersect(rownames(cortex_mtx),rownames(gut_rel_mtx))
#merge the matrices, change colnames to include annotation and write the files for cell signal analysis
merged_mtx=cbind(gut_rel_mtx[comm_genes,], cortex_mtx[comm_genes,])
colnames(merged_mtx)=c(paste0(gut_srat@meta.data$final_annot,":",colnames(gut_rel_mtx)),
                       paste0("Cortex:", colnames(cortex_mtx)))

writeMM(merged_mtx, "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_gut/sc_gut_with_cortex/gut.mtx")
write.table(rownames(merged_mtx), "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_gut/sc_gut_with_cortex/gut_rowNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(colnames(merged_mtx), "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_gut/sc_gut_with_cortex/gut_columnNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)

#########set up bulk files##############
#TCGA colon and rectum paths
rc2 = readRDS('/lustre/scratch119/realdata/mdt1/team274/my4/bulkRNAseq/rse_gene_TCGA.RDS')
len=read.table('matts-bits-of-code-master/cellSignalAnalysisV2/bulkData/SangerProjectTARGET_ALL_phase1/SangerProjectTARGET_ALL_phase1_SRR452447.tsv', sep = "\t", header = T)
obs=read.csv("onecs/obs.csv")
mDat = colData(rc2)

mDat$UniqueSampleID = rownames(mDat)                                                                                                                                                                                                                                                                                                                                                                                                                                                     
mDat$age = mDat$cgc_case_age_at_diagnosis*12
mDat$biopsy = mDat$gdc_cases.samples.sample_type
mDat$tissue = as.character(mDat$xml_tumor_tissue_site)

mDat_col=mDat[which(mDat$tissue%in%c("Rectum", "Colon")),]
mDat_col=mDat_col[,colnames(mDat_col)[which(!colnames(mDat_col)%in%names(which(apply(mDat_col, 2, function(x){sum(is.na(x))})==715)))]]
saveRDS(mDat_col,'mDat_gut.rds')

paths=paste0("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/CellDeconvolution/Results/bulkRNAseq/TCGA/frag/TCGA_", rownames(mDat_col), ".tsv")

write.table(paths, "colon_bulk_paths.txt", sep = "\t", row.names = F, col.names = F, quote = F)

#fetal bulk files 

load('Foetal_RawCounts_Zilbauer.RData')

zilb_counts=counts
zilb_meta=foetal.info

zilb_meta$sample.id[which(zilb_meta$sample.type=="purified")]

zilb_fetal_counts=zilb_counts[,zilb_meta$sample.id[which(zilb_meta$sample.type=="purified")]]

exmp=read.table(
  'matts-bits-of-code-master/cellSignalAnalysis/bulkData/SangerProjectTARGET_ALL_phase1/SangerProjectTARGET_ALL_phase1_SRR452447.tsv',
  sep = "\t", header = T)
rownames(exmp)=exmp$geneName
exp1=intersect(rownames(exmp), rownames(zilb_counts))
exmp=exmp[exp1,]
zilb_fetal_counts=zilb_fetal_counts[exp1,]

zilb_fetal_counts$geneLengths=exmp$geneLengths
zilb_fetal_counts$geneName=rownames(zilb_fetal_counts)

for (x in zilb_meta$sample.id[which(zilb_meta$sample.type=="purified")]) {
  a = zilb_fetal_counts[,c("geneName","geneLengths",x)]
  write.table(a, paste0('cellsig_gut/fetal_bulk/',x,'.tsv'), sep="\t", row.names = F, col.names = T, quote = F)
}

# gtex brain as a negative control - change files to set gene lengths to 1 for everything
gtex_meta=read.delim('/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/CellDeconvolution/Results/bulkRNAseq/GTEX/Metadata.tsv')

rel_ids=gtex_meta$UniqueSampleID[which(gtex_meta$Tissue%in%c("Brain - Frontal Cortex (BA9)",
                                                             "Brain - Cerebellar Hemisphere","Brain - Cortex"))]
rel_gtex_meta=gtex_meta[which(gtex_meta$Tissue%in%c("Brain - Frontal Cortex (BA9)",
                                                    "Brain - Cerebellar Hemisphere","Brain - Cortex")),]
gtex_paths=paste0('/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/CellDeconvolution/Results/bulkRNAseq/GTEX/frag/',rel_ids, '.tsv')
gtex_tissues=data.frame(path=gtex_paths, tissue=rel_gtex_meta$Tissue)
# example of how files need to be structured
bla=read.table("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/CellDeconvolution/Results/bulkRNAseq/GTEX/frag/GTEX_SRR627421.tsv", sep = "\t", header=T)

#set gtex gene lengths to 1 

for (i in gtex_paths) {
  file=read.table(i, sep = "\t", header=T)
  file$geneLengths=1
  
  write.table(file, paste0('cellsig_gut/gtex_brain_bulk/', basename(i)), sep = "\t", col.names = T, row.names = F, quote = F)
  
}

new_paths=paste0('/home/jovyan/Dediff/cellsig_gut/gtex_brain_bulk/', basename(gtex_paths))
write.table(new_paths, 'cellsig_gut/gtex_brain_paths.txt', sep = "\t", quote=F, row.names = F, col.names = F)
write.table(gtex_tissues, 'gtex_brain_tissues.txt', sep = "\t", row.names = F, col.names = T, quote=F)

#boardman colon data

boardman_meta=read.table('boardman_data/boardman_meta.txt', sep = ",", header = T)
boardman_data=read.delim('boardman_data/boardman_data.tsv', sep = "\t", header = T)

boardman_rel=boardman_data[,c(4,5,grep("count", colnames(boardman_data)))]

colnames(boardman_rel)=gsub("_count", "", colnames(boardman_rel))

length(which(boardman_meta$sample%in%colnames(boardman_rel)))

colnames(boardman_rel)=gsub("Unk.", "", colnames(boardman_rel))

boardman_meta_new=data.frame(sample=colnames(boardman_rel)[3:81])
which(boardman_meta_new$sample%in%boardman_meta$sample)
boardman_meta_new$sample_type="x"

boardman_meta_new$sample_type=boardman_meta$sample_type[match(boardman_meta_new$sample, boardman_meta$sample)]

boardman_meta_new$sample_type[which(is.na(boardman_meta_new$sample_type))]="unknown"

colnames(boardman_rel)[1]="geneLengths"
colnames(boardman_rel)[2]="geneName"

for (x in colnames(boardman_rel)[3:81]) {
  write.table(boardman_rel[1:64253,c("geneName","geneLengths",x)], paste0('cellsig_gut/boardman_data/',x,'.tsv'), sep="\t", row.names = F, col.names = T, quote = F)
}

boardman_paths=paste0('/home/jovyan/Dediff/cellsig_gut/boardman_data/', colnames(boardman_rel)[3:81], '.tsv')

write.table(boardman_paths, "/home/jovyan/Dediff/cellsig_gut/boardman_paths.txt",sep = "\t", row.names = F, col.names = F, quote=F)

# at this point, select relevant paths and run cell signal analysis [Matt's code]



