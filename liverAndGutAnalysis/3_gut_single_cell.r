source('process_seurat.R')
source('logisticRegression.R')
library(stringr)
library(circlize)
gut_srat=readRDS('/lustre/scratch117/casm/team274/gk14/rasa_final_object_soupx/rasa_rel_srat.rds')

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
table(gut_srat$final_annot)
gut_srat=subset(gut_srat, final_annot%in%c("Fetal_CLDN10+ cells","Fetal_Microfold cell", "Fetal_Paneth", "Fetal_Tuft"), invert=T)
#subset cortex from adrenal gland dataset
adr_all=readRDS('/lustre/scratch117/casm/team274/gk14/cellxgene/adr_all.rds')
cortex=subset(adr_all, Annotation%in%c("Cortex"))
cortex@meta.data$final_annot=cortex@meta.data$Annotation
#relabel the gene names to match
gut_genes=read.csv('onecs/obs.csv')
adr_genes=read.table("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/SCP/Data/fAdrenal19wk/cellranger302_count_32644_WSSS_F_Adr8710632_GRCh38-1_2_0/filtered_feature_bc_matrix/features.tsv.gz", sep = "\t", header = F)

gut_mtx=gut_srat@assays$RNA@counts
cortex_mtx=cortex@assays$RNA@counts

rownames(gut_mtx)=gut_genes$gene_ids
rownames(cortex_mtx)=adr_genes$V1

inter=intersect(rownames(gut_mtx), rownames(cortex_mtx))
rownames(gut_genes)=gut_genes$gene_ids

inter_with_ee_gast=intersect(intersect(intersect(rownames(gut_mtx),
                                                 rownames(cortex_mtx)),
                                       rownames(early_emb_mtx)),
                             rownames(gastrulation_mtx))

merged_plus_early_stuff=cbind(gut_mtx[inter_with_ee_gast,], cortex_mtx[inter_with_ee_gast,],
                              early_emb_mtx[inter_with_ee_gast,],gastrulation_rpk[inter_with_ee_gast,])
rownames(merged_plus_early_stuff)=adr_genes[rownames(merged_plus_early_stuff),]$V2

gut_ee_gast_colnames=c(gut_srat@meta.data$final_annot, cortex@meta.data$final_annot, as.character(Idents(early_emb)), gast_meta_germl$wide_clust)

gut_gast_ee_ref=CreateSeuratObject(merged_plus_early_stuff)
gut_gast_ee_ref@meta.data$annot=gut_ee_gast_colnames
saveRDS(gut_gast_ee_ref, '/lustre/scratch117/casm/team274/gk14/Dediff/gut_gast_ee_ref.rds')
gut_gast_ee=readRDS('/lustre/scratch117/casm/team274/gk14/Dediff/gut_gast_ee_ref.rds')
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



#lee colorec data

  
crc_lee=read.table('GSE132465/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt', header = T, sep = "\t", row.names = 1)
crc_lee_annot=read.table('GSE132465/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt', header = T, sep = "\t")

crc_lee=as.matrix(crc_lee)
crc_lee=Matrix(crc_lee, sparse=T)
crc_lee_annot$Index=str_replace(crc_lee_annot$Index, pattern = "-", replacement = ".")
crc_lee_annot_tum=crc_lee_annot[which(crc_lee_annot$Class=="Tumor"),]
dim(crc_lee_annot)


crc_lee_tum=crc_lee[,crc_lee_annot_tum$Index]

rownames(crc_lee_annot_tum)=crc_lee_annot_tum$Index

crc_lee_srat=CreateSeuratObject(crc_lee_tum, meta.data = crc_lee_annot_tum)

crc_lee_srat=process_seurat(crc_lee_srat)
crc_lee_srat@meta.data$Cell_subtype
crc_lee_srat@meta.data$wide_ann_pt=crc_lee_srat@meta.data$Cell_type
crc_lee_srat@meta.data$wide_ann_pt[which(crc_lee_srat@meta.data$wide_ann_pt=="Epithelial cells")]=paste0("Epithelial_cells_", crc_lee_srat@meta.data$Patient[which(crc_lee_srat@meta.data$wide_ann_pt=="Epithelial cells")])

crc_lee_srat@meta.data$wide_ann_pt[which(crc_lee_srat@meta.data$wide_ann_pt=="Epithelial cells")]=paste0(crc_lee_srat@meta.data$Cell_subtype[which(crc_lee_srat@meta.data$wide_ann_pt=="Epithelial cells")], crc_lee_srat@meta.data$Patient[which(crc_lee_srat@meta.data$wide_ann_pt=="Epithelial cells")])

saveRDS(crc_lee_srat,'/lustre/scratch117/casm/team274/gk14/Dediff/GSE132465/crc_lee_srat.rds')
crc_lee_srat=readRDS('/lustre/scratch117/casm/team274/gk14/Dediff/GSE132465/crc_lee_srat.rds')

tum_inter=intersect(rownames(new_g_merged), rownames(crc_lee_srat))
tum_ee_inter=intersect(rownames(gut_gastr_ee), rownames(crc_lee_srat))

tr_gut1=trainModel(new_g_merged@assays$RNA@counts[tum_inter, ], new_g_merged@meta.data$final_annot, workers = NULL, minCells = 450)
tr_gut2=trainModel(new_g_merged@assays$RNA@counts[tum_inter, ], new_g_merged@meta.data$final_annot, workers = NULL, minCells = 450)
tr_gut3=trainModel(new_g_merged@assays$RNA@counts[tum_inter, ], new_g_merged@meta.data$final_annot, workers = NULL, minCells = 450)

tr_gut_fullref=trainModel(gut_gastr_ee@assays$RNA@counts[tum_ee_inter,], gut_gastr_ee@meta.data$annot, workers = NUL, minCells = 450)

pr_gut1=predictSimilarity(tr_gut1, crc_lee_srat@assays$RNA@counts[tum_inter, ], crc_lee_srat@meta.data$wide_ann_pt, logits = F)
similarityHeatmap(pr_gut1)

pr_fullref=predictSimilarity(tr_gut_fullref, crc_lee_srat@assays$RNA@counts[tum_ee_inter, ], crc_lee_srat@meta.data$wide_ann_pt, logits = F)
similarityHeatmap(pr_gut1)
similarityHeatmap(pr_fullref)

cell_counts=data.frame(unclass(table(crc_lee_srat@meta.data$Cell_subtype, crc_lee_srat@meta.data$Patient)))
cms_per_patient_count=cell_counts[c("CMS1","CMS2","CMS3","CMS4"),]

rownames(cms_per_patient_count)[apply(cms_per_patient_count,2,which.max)]

split_dat=data.frame(rname=rownames(pr_gut1),split=c("B cells", rownames(cms_per_patient_count)[apply(cms_per_patient_count,2,which.max)],
                                                     "Mast cells", "Myeloids","Stromal cells", "T cells" ))


similarityHeatmap(pr_gut1, row_split=split_dat$split)

col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))

Heatmap(pr_gut1, cluster_columns=F, cluster_rows = F, col=col_fun,  column_gap = unit(3, "mm"), row_names_side = "left",
        row_split=split_dat$split,
        width = ncol(pr_gut1)*unit(4, "mm"), 
        height = nrow(pr_gut1)*unit(4, "mm"))
pr_gut2=predictSimilarity(tr_gut1, crc_lee_srat@assays$RNA@counts[tum_inter, ], crc_lee_srat@meta.data$Cell_subtype, logits = F)
Heatmap(pr_gut2, cluster_columns=F, cluster_rows = F, col=col_fun,  column_gap = unit(3, "mm"), row_names_side = "left",
        width = ncol(pr_gut1)*unit(4, "mm"), 
        height = nrow(pr_gut1)*unit(4, "mm"))


Heatmap(pr_fullref, cluster_columns=F, cluster_rows = F, col=col_fun,  column_gap = unit(3, "mm"), row_names_side = "left",
        width = ncol(pr_fullref)*unit(4, "mm"), 
        height = nrow(pr_fullref)*unit(4, "mm"))

DimPlot(crc_lee_srat, group.by = "wide_ann_pt")
crc_lee_srat@meta.data$wide_ann_pt2=crc_lee_srat@meta.data$wide_ann_pt
crc_lee_srat@meta.data$wide_ann_pt2[grep("Epithelial",crc_lee_srat@meta.data$wide_ann_pt2)]="Tumor"
colfunc <- colorRampPalette(c("#084594", "#C6DBEF"))
tol9qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")

DimPlot(crc_lee_srat, group.by = "wide_ann_pt", label=T,
        cols = c( "#771122",colfunc(23),"#771122","#DDCC77","#44AA99","#7E481C"),
        pt.size = 0.1) +theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
                              axis.ticks = element_blank(), legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1") +ggtitle("")
pdf("crc_umap.pdf", height = 3, width=3,useDingbats = F)
gg=DimPlot(crc_lee_srat, group.by = "wide_ann_pt", label=T,
        cols = c( "#771122",colfunc(23),"#771122","#DDCC77","#44AA99","#7E481C"),
        pt.size = 0.1) +theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
                              axis.ticks = element_blank(), legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1") +ggtitle("")
plot(gg)
dev.off()


gut_colsplit_df=data.frame(celltype=colnames(pr_fullref))
gut_colsplit_df$wide_annot=gut_colsplit_df$celltype
gut_colsplit_df$wide_annot[which(gut_colsplit_df$wide_annot%in%c("Hypoblast", "Epiblast"))]="Early embryo"
gut_colsplit_df$wide_annot[which(gut_colsplit_df$wide_annot%in%c("Mesoderm",
                                                                 "Endoderm",
                                                                 "Non-Neural Ectoderm"))]="Gastrulation"


r_ct=read.csv('/lustre/scratch117/casm/team274/gk14/rasa_final_object_soupx/var.csv')

gut_colsplit_df$wide_annot[which(gut_colsplit_df$wide_annot%in%c("Fetal_BEST4+ epithelial",
                                                                 "Fetal_Colonocyte","Fetal_Distal progenitor", 
                                                                 "Fetal_Enterocyte",
                                                                 "Fetal_Goblet cell", 
                                                                 "Fetal_Proximal progenitor",
                                                                 "Fetal_Stem cells" ,
                                                                 "Fetal_TA"))]="Fetal_epithelial"

gut_colsplit_df$wide_annot[which(gut_colsplit_df$wide_annot%in%c("Adult_Paneth",
                                                                 "Adult_BEST2+ Goblet cell",
                                                                 "Adult_BEST4+ epithelial",
                                                                 "Adult_Colonocyte", 
                                                                 "Adult_Goblet cell",
                                                                 "Adult_Microfold cell",
                                                                 "Adult_Stem cells",
                                                                 "Adult_Tuft",
                                                                 "Adult_TA"))]="Adult_epithelial"



gut_colsplit_df$wide_annot[which(gut_colsplit_df$wide_annot%in%c("Adult_B cells",
                                                                 "Adult_Myeloid",
                                                                 "Adult_Plasma cells",
                                                                 "Adult_T cells"))]="Adult_blood/immune"


gut_colsplit_df$wide_annot[which(gut_colsplit_df$wide_annot%in%c("Fetal_B cells" , "Fetal_Red blood cells",
                                                                 "Fetal_T cells", "Fetal_Myeloid" ))]="Fetal_blood/immune"

gut_colsplit_df$wide_annot[which(gut_colsplit_df$wide_annot%in%c("Adult_Enteroendocrine",
                                                                 "Adult_Endothelial",
                                                                 "Adult_Mesenchymal",
                                                                 "Adult_Neuronal" ))]="Adult_Other"




gut_colsplit_df$wide_annot[which(gut_colsplit_df$wide_annot%in%c("Fetal_Endothelial",
                                                                 "Fetal_Enteroendocrine",
                                                                 "Fetal_Neuronal",
                                                                 "Fetal_Mesenchymal"))]="Fetal_Other"


gut_colsplit_df$wide_annot=factor(gut_colsplit_df$wide_annot, levels=c("Early embryo" , 
                                                                       "Gastrulation" ,
                                                                       "Adult_epithelial" ,
                                                                       "Fetal_epithelial",
                                                                       "Adult_blood/immune",
                                                                       "Fetal_blood/immune",
                                                                       "Adult_Other",
                                                                       "Fetal_Other",
                                                                       "Cortex"))

gut_row_ord=rownames(pr_fullref)[c(2:28,1)]

gut_colsplit_df$celltype=factor(gut_colsplit_df$celltype, levels=c(
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
  "Cortex"))

tr_gut_fullref=trainModel(gut_gastr_ee@assays$RNA@counts[tum_ee_inter,], gut_gastr_ee@meta.data$annot, workers = NULL, minCells = 300)
pr_fullref=predictSimilarity(tr_gut_fullref, crc_lee_srat@assays$RNA@counts[tum_ee_inter, ], crc_lee_srat@meta.data$wide_ann_pt, logits = F)

Heatmap(pr_fullref[,as.character(gut_colsplit_df$celltype)], cluster_columns=F, cluster_rows = F, col=col_fun, row_order=gut_row_ord,
        column_split=gut_colsplit_df$wide_annot, column_order=levels(gut_colsplit_df$celltype),
        column_gap = unit(3, "mm"), row_names_side = "left",
        width = ncol(pr_fullref)*unit(4, "mm"), 
        height = nrow(pr_fullref)*unit(4, "mm"))


tr_gut_fullref2=trainModel(gut_gastr_ee@assays$RNA@counts[tum_ee_inter,], gut_gastr_ee@meta.data$annot, workers = NULL, minCells = 300)
pr_fullref2=predictSimilarity(tr_gut_fullref2, crc_lee_srat@assays$RNA@counts[tum_ee_inter, ], crc_lee_srat@meta.data$wide_ann_pt, logits = F)

Heatmap(pr_fullref2[,as.character(gut_colsplit_df$celltype)], cluster_columns=F, cluster_rows = F, col=col_fun, row_order=gut_row_ord,
        column_split=gut_colsplit_df$wide_annot, column_order=levels(gut_colsplit_df$celltype),
        column_gap = unit(3, "mm"), row_names_side = "left",
        width = ncol(pr_fullref)*unit(4, "mm"), 
        height = nrow(pr_fullref)*unit(4, "mm"))

tr_gut_fullref3=trainModel(gut_gastr_ee@assays$RNA@counts[tum_ee_inter,], gut_gastr_ee@meta.data$annot, workers = NULL, minCells = 300)
pr_fullref3=predictSimilarity(tr_gut_fullref3, crc_lee_srat@assays$RNA@counts[tum_ee_inter, ], crc_lee_srat@meta.data$wide_ann_pt, logits = F)

Heatmap(pr_fullref3[,as.character(gut_colsplit_df$celltype)], cluster_columns=F, cluster_rows = F, col=col_fun, row_order=gut_row_ord,
        column_split=gut_colsplit_df$wide_annot, column_order=levels(gut_colsplit_df$celltype),
        column_gap = unit(3, "mm"), row_names_side = "left",
        width = ncol(pr_fullref)*unit(4, "mm"), 
        height = nrow(pr_fullref)*unit(4, "mm"))



