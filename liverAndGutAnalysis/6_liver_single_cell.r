library(SoupX)
library(Seurat)
source('logisticRegression.R')
source('process_seurat.R')
#####liver single cell validation######
merged_liver_srat=readRDS('/lustre/scratch117/casm/team274/gk14/Dediff/merged_liver/merged_liver_srat.rds')
adr_all=readRDS('/lustre/scratch117/casm/team274/gk14/cellxgene/adr_all.rds')
cortex=subset(adr_all, Annotation%in%c("Cortex"))
cortex_new=cortex@assays$RNA@counts[rownames(merged_liver_srat),]
cortex_new=CreateSeuratObject(cortex_new, meta.data = cortex@meta.data)
cortex_new@meta.data$annot=cortex_new@meta.data$Annotation
liver_and_cortex=merge(merged_liver_srat, cortex_new)

liver_and_cortex = NormalizeData(liver_and_cortex)
liver_and_cortex =FindVariableFeatures(liver_and_cortex, selection.method = "vst", nfeatures = 2000)
liver_and_cortex = CellCycleScoring(liver_and_cortex, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

liver_and_cortex@meta.data$wide_annot=liver_and_cortex@meta.data$annot
liver_and_cortex@meta.data$wide_annot[which(liver_and_cortex@meta.data$wide_annot%in%c("Fetal_Mono-Mac","Fetal_B-cells","Fetal_DCs",
                                                                                       "Fetal_Monocyte",
                                                                                       "Fetal_Early lymphoid_T lymphocyte",
                                                                                       "Fetal_VCAM1+ EI macrophage", "Fetal_ILC precursor",
                                                                                       "Fetal_Neutrophil-myeloid progenitor",
                                                                                       "Fetal_Monocyte precursor", 
                                                                                       "Fetal_NK",
                                                                                       "Fetal_Mast cell", 
                                                                                       "Fetal_Kupffer Cell"))]="Fetal_leukocytes"
liver_and_cortex@meta.data$wide_annot[which(liver_and_cortex@meta.data$wide_annot%in%c("Adult_Non-inflammatory_macs", 
                                                                                       "Adult_T-cells", 
                                                                                       "Adult_Inflammatory_macs",
                                                                                       "Adult_NK-like_cells",
                                                                                       "Adult_Mature_Bcells"))]="Adult_leukocytes"
liver_and_cortex@meta.data$wide_annot[which(liver_and_cortex@meta.data$wide_annot%in%c("Adult_Hepatocytes", 
                                                                                       "Fetal_Hepatocyte", 
                                                                                       "Adult_HHyP",
                                                                                       "Fetal_HHyP",
                                                                                       "Adult_BECs"))]="EPCAM+"

liver_and_cortex@meta.data$wide_annot[which(liver_and_cortex@meta.data$wide_annot%in%c("Fetal_Erythroid", 
                                                                                       "Fetal_MEMP", 
                                                                                       "Fetal_Megakaryocyte",
                                                                                       "Fetal_MEMP",
                                                                                       "Fetal_HSC_MPP",
                                                                                       "Adult_Plasma_cells",
                                                                                       "Adult_Erythroid_cells"))]="Other blood cells"
liver_and_cortex@meta.data$wide_annot[which(liver_and_cortex@meta.data$wide_annot%in%c("Cortex"))]="Negative control"

liver_and_cortex@meta.data$wide_annot[which(liver_and_cortex@meta.data$wide_annot%in%c("Fetal_Endothelial cell",
                                                                                       "Adult_LSECs",
                                                                                       "Adult_Portal_endothelium"))]="Endothelium"
liver_and_cortex@meta.data$wide_annot[which(liver_and_cortex@meta.data$wide_annot%in%c("Fetal_leukocytes" ,
                                                                                       "Adult_leukocytes"))]="Leukocytes"

saveRDS(liver_and_cortex, '/lustre/scratch117/casm/team274/gk14/Dediff/liver_and_cortex.rds')

merged_liver_srat=readRDS('/lustre/scratch117/casm/team274/gk14/Dediff/merged_liver/merged_liver_srat.rds')

length(intersect(rownames(merged_liver_srat), ensembl_adrenal$V2))

l_and_c=cbind(merged_liver_srat@assays$RNA@counts, cortex@assays$RNA@counts[rownames(merged_liver_srat),])
ens_liver=ensembl_adrenal[rownames(merged_liver_srat),]

colnames(l_and_c)=c(paste0(merged_liver_srat@meta.data$annot,":",colnames(merged_liver_srat)),
                    paste0("Cortex:", colnames(cortex@assays$RNA@counts)))
celltypes_liver=c(unique(merged_liver_srat@meta.data$annot), "Cortex")

rownames(l_and_c)=ens_liver$V1
inter_liver_with_ee_gast=intersect(intersect(rownames(l_and_c),
                                            rownames(early_emb_mtx)),
                                  rownames(gastrulation_mtx))

liver_gastr_ee=cbind(l_and_c[inter_liver_with_ee_gast,],
                    early_emb_mtx[inter_liver_with_ee_gast,],gastrulation_rpk[inter_liver_with_ee_gast,])
rownames(liver_gastr_ee)=adr_genes[rownames(liver_gastr_ee),]$V2
liver_gastr_ee_annot=c(merged_liver_srat@meta.data$annot, cortex@meta.data$Annotation, as.character(Idents(early_emb)), gast_meta_germl$wide_clust)

liver_gastr_ee=CreateSeuratObject(liver_gastr_ee)
liver_gastr_ee@meta.data$annot=liver_gastr_ee_annot
saveRDS(liver_gastr_ee, '/lustre/scratch117/casm/team274/gk14/Dediff/liver_gast_ee_ref.rds')
liver_gastr_ee=readRDS('/lustre/scratch117/casm/team274/gk14/Dediff/liver_gast_ee_ref.rds')

#Get the inhouse HB data

paths=c('/lustre/scratch117/casm/team274/my4/oldScratch/Data/Disease/Hepatoblastoma/GOSH020/cellranger302_count_27838_4602STDY7741866_GRCh38-1_2_0/filtered_feature_bc_matrix/',
        '/lustre/scratch117/casm/team274/my4/oldScratch/Data/Disease/Hepatoblastoma/GOSH020/cellranger302_count_27838_4602STDY7741867_GRCh38-1_2_0/filtered_feature_bc_matrix/',
        '/lustre/scratch117/casm/team274/my4/oldScratch/Data/Disease/Hepatoblastoma/GOSH020/cellranger302_count_27838_4602STDY7741868_GRCh38-1_2_0/filtered_feature_bc_matrix/',
        '/lustre/scratch117/casm/team274/my4/oldScratch/Data/Disease/Hepatoblastoma/GOSH020/cellranger302_count_27926_4602STDY7762295_GRCh38-1_2_0/filtered_feature_bc_matrix/',
        '/lustre/scratch117/casm/team274/my4/oldScratch/Data/Disease/Hepatoblastoma/GOSH020/cellranger302_count_27926_4602STDY7762297_GRCh38-1_2_0/filtered_feature_bc_matrix/',
        '/lustre/scratch117/casm/team274/my4/oldScratch/Data/Disease/Hepatoblastoma/GOSH020/cellranger302_count_27926_4602STDY7762299_GRCh38-1_2_0/filtered_feature_bc_matrix/')

inhouse_hb_mtx=Read10X(paths)

inhouse_hb_srat=CreateSeuratObject(inhouse_hb_mtx)
inhouse_hb_srat[["percent.mt"]] = PercentageFeatureSet(inhouse_hb_srat, pattern = "^MT-")

VlnPlot(inhouse_hb_srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

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

DimPlot(inhouse_hb_srat, label=T)


FeaturePlot(inhouse_hb_srat, "percent.mt")


lit_genes=c('PTPRC','PDGFRB','EPCAM','PECAM1','HBA1','HBB')

FeaturePlot(inhouse_hb_srat,lit_genes)
DimPlot(inhouse_hb_srat, label=T)
labAnch_inh = FindTransferAnchors(reference = liver_and_cortex,query = inhouse_hb_srat,dims=seq(30))
pred_inh = TransferData(labAnch_inh,refdata=liver_and_cortex@meta.data$annot,dims=seq(30))
inhouse_hb_srat = AddMetaData(inhouse_hb_srat,metadata=pred_inh)
DimPlot(inhouse_hb_srat, group.by = "predicted.id", label = T)

marks=quickMarkers(inhouse_hb_srat@assays$RNA@counts, inhouse_hb_srat@meta.data$seurat_clusters, N=10)

FeaturePlot(inhouse_hb_srat,c('prediction.score.Cortex', 'prediction.score.Adult_Erythroid_cells'))
DimPlot(inhouse_hb_srat, label=T)

#annotate the cells based on marker expression and automated annot

inhouse_hb_srat@meta.data$annot="x"
inhouse_hb_srat@meta.data$annot[which(inhouse_hb_srat@meta.data$seurat_clusters%in%c(8,9,10,11,13,14,18,19,20,21,22))]="Leukocytes"
inhouse_hb_srat@meta.data$annot[which(inhouse_hb_srat@meta.data$seurat_clusters%in%c(12))]="Erythroid cells"
inhouse_hb_srat@meta.data$annot[which(inhouse_hb_srat@meta.data$seurat_clusters%in%c(4))]="LSECs"
inhouse_hb_srat@meta.data$annot[which(inhouse_hb_srat@meta.data$seurat_clusters%in%c(15))]="Hepatic stellate cells"
inhouse_hb_srat@meta.data$annot[which(inhouse_hb_srat@meta.data$seurat_clusters%in%c(0,1,2,3,5,6,7,16,17))]="Cancer cells"
inhouse_hb_srat@meta.data$annot[which(inhouse_hb_srat@meta.data$seurat_clusters%in%c(23))]="Unknown"
DimPlot(inhouse_hb_srat, group.by = "annot")

inhouse_hb_srat@meta.data$annot=factor(inhouse_hb_srat@meta.data$annot,
                                       levels=c("Cancer cells","Hepatic stellate cells", 
                                                "LSECs", "Erythroid cells",
                                                "Leukocytes","Unknown"))

colfunc <- colorRampPalette(c("#084594", "#C6DBEF")) 

pdf("hb_inhouse_umap.pdf", height = 3, width=3,useDingbats = F)
gg=DimPlot(inhouse_hb_srat, group.by = "annot", label=T,
           cols = c("#084594","#7E481C","#CC6677" , "#771122","#DDCC77", "lightgrey"),
           pt.size = 0.1) +theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
                                 axis.ticks = element_blank(), legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1") +ggtitle("")
plot(gg)
dev.off()


inhouse_shared_genes=intersect(rownames(inhouse_hb_srat), rownames(liver_and_cortex))

inhouse_ee_gast_shared_genes=intersect(rownames(inhouse_hb_srat), rownames(liver_gastr_ee))
source('logisticRegression.R')
tr_inh=trainModel(liver_and_cortex@assays$RNA@counts[inhouse_shared_genes, ],
                  as.character(liver_and_cortex@meta.data$annot),
                  workers=NULL, minCells = 300)

tr_inh_nomin=trainModel(liver_and_cortex@assays$RNA@counts[inhouse_shared_genes, ],
                  as.character(liver_and_cortex@meta.data$annot),
                  workers=NULL)
pr_inh=predictSimilarity(tr_inh_nomin,
                         inhouse_hb_srat@assays$RNA@counts[inhouse_shared_genes,], logits = T)


tr_inh_ee_gast=trainModel(liver_gastr_ee@assays$RNA@counts[inhouse_ee_gast_shared_genes, ],
                  as.character(liver_gastr_ee@meta.data$annot),
                  workers=NULL, minCells = 300)
pr_inh_ee_gast=predictSimilarity(tr_inh_ee_gast,
                                 inhouse_hb_srat@assays$RNA@counts[inhouse_ee_gast_shared_genes,], 
                                 inhouse_hb_srat@meta.data$annot,logits = T)

similarityHeatmap(pr_inh_ee_gast)

Heatmap(pr_inh, cluster_columns=F, cluster_rows = F, column_split=colsplit_df$wide_annot, col=col_fun, 
        column_order=levels(colsplit_df$celltype),  column_gap = unit(3, "mm"), row_names_side = "left",
        column_title_side = "bottom", width = ncol(pr_inh)*unit(4, "mm"), 
        height = nrow(pr_inh)*unit(4, "mm"))
unique(colsplit_df$wide_annot)
colsplit_df=data.frame(celltype=colnames(pr_inh))
colsplit_df$wide_annot=colsplit_df$celltype
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%c("Fetal_Mono-Mac","Fetal_B-cells","Fetal_DCs",
                                                                                       "Fetal_Monocyte",
                                                                                       "Fetal_Early lymphoid_T lymphocyte",
                                                                                       "Fetal_VCAM1+ EI macrophage", "Fetal_ILC precursor",
                                                                                       "Fetal_Neutrophil-myeloid progenitor",
                                                                                       "Fetal_Monocyte precursor", 
                                                                                       "Fetal_NK",
                                                                                       "Fetal_Mast cell", 
                                                                                       "Fetal_Kupffer Cell"))]="Fetal_leukocytes"
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%c("Adult_Non-inflammatory_macs", 
                                                                                       "Adult_T-cells", 
                                                                                       "Adult_Inflammatory_macs",
                                                                                       "Adult_NK-like_cells",
                                                                                       "Adult_Mature_Bcells"))]="Adult_leukocytes"
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%c("Adult_Hepatocytes", 
                                                                                       "Fetal_Hepatocyte", 
                                                                                       "Adult_HHyP",
                                                                                       "Fetal_HHyP",
                                                                                       "Adult_BECs"))]="EPCAM+"

colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%c("Fetal_Erythroid", 
                                                                                       "Fetal_MEMP", 
                                                                                       "Fetal_Megakaryocyte",
                                                                                       "Fetal_MEMP",
                                                                                       "Fetal_HSC_MPP",
                                                                                       "Adult_Plasma_cells",
                                                                                       "Adult_Erythroid_cells"))]="Other blood cells"
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%c("Cortex"))]="Negative control"
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%c("Fetal_Endothelial cell",
                                                                                       "Adult_LSECs",
                                                                                       "Adult_Portal_endothelium"))]="Endothelium"

colsplit_df$wide_annot=factor(colsplit_df$wide_annot, levels=c("EPCAM+", "Adult_Hepatic_stellate_cells", "Endothelium",
                                                               "Other blood cells","Adult_leukocytes", "Fetal_leukocytes", 
                                                               "Fetal_Fibroblast",  "Negative control" ))
colsplit_df$celltype=factor(colsplit_df$celltype, levels = c("Fetal_Hepatocyte",
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
                                                             "Cortex"))
tr_inh_wide_annot=trainModel(liver_and_cortex@assays$RNA@counts[inhouse_shared_genes, ],
                             as.character(liver_and_cortex@meta.data$wide_annot), workers=NULL)
pr_inh_wide_annot=predictSimilarity(tr_inh_wide_annot,inhouse_hb_srat@assays$RNA@counts[inhouse_shared_genes,])
col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
similarityHeatmap(pr_inh_wide_annot, cluster_rows=T, use_raster=F, row_split=inhouse_hb_srat@meta.data$annot)




theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
      axis.ticks = element_blank(), legend.position = "none")

hcc_sc_paths= c('/lustre/scratch117/casm/team274/gk14/Dediff/ho_hcc_sc/SRR14424777/output/Gene/filtered/', 
                '/lustre/scratch117/casm/team274/gk14/Dediff/ho_hcc_sc/SRR14424778/output/Gene/filtered/', 
                '/lustre/scratch117/casm/team274/gk14/Dediff/ho_hcc_sc/SRR14424779/output/Gene/filtered/',
                '/lustre/scratch117/casm/team274/gk14/Dediff/ho_hcc_sc/SRR14424780/output/Gene/filtered/',
                '/lustre/scratch117/casm/team274/gk14/Dediff/ho_hcc_sc/SRR14424781/output/Gene/filtered/', 
                '/lustre/scratch117/casm/team274/gk14/Dediff/ho_hcc_sc/SRR14424782/output/Gene/filtered/',
                '/lustre/scratch117/casm/team274/gk14/Dediff/ho_hcc_sc/SRR14424783/output/Gene/filtered/',
                '/lustre/scratch117/casm/team274/gk14/Dediff/ho_hcc_sc/SRR14424784/output/Gene/filtered/')
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
DimPlot(hcc_sc_srat, label=T)

hcc_inter=intersect(rownames(hcc_sc_srat), rownames(liver_and_cortex))
train_for_hcc=trainModel(liver_and_cortex@assays$RNA@counts[hcc_inter,], as.character(liver_and_cortex@meta.data$annot), workers=NULL)
hcc_ps=predictSimilarity(train_for_hcc, hcc_sc_srat@assays$RNA@counts[hcc_inter,], hcc_sc_srat@meta.data$annot2, logits = F)

hcc_ps_orig_ident=predictSimilarity(train_for_hcc, hcc_sc_srat@assays$RNA@counts[hcc_inter,], hcc_sc_srat@meta.data$annot3, logits = F)
similarityHeatmap(hcc_ps)
Heatmap(hcc_ps, column_order=levels(colsplit_df$celltype))
row_ord=c("cancer_1","cancer_2","cancer_23","cancer_28","cancer_31","cancer_32","cancer_33",
           "cancer_15","cancer_17","cancer_29","cancer_35",
          "cancer_19","cancer_7","cancer_6","cancer_0","cancer_22",
          "cancer_26", "Endothelium","Leukocytes")
col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
Heatmap(hcc_ps, cluster_columns=F, cluster_rows = F, column_split=colsplit_df$wide_annot, col=col_fun, 
        column_order=levels(colsplit_df$celltype), row_order=row_ord, column_gap = unit(3, "mm"), row_names_side = "left",
        column_title_side = "bottom", width = ncol(hcc_ps)*unit(4, "mm"), 
        height = nrow(hcc_ps)*unit(4, "mm"))

Heatmap(hcc_ps_orig_ident, cluster_columns=F, cluster_rows = F, column_split=colsplit_df$wide_annot, col=col_fun, 
        column_order=levels(colsplit_df$celltype), column_gap = unit(3, "mm"), row_names_side = "left",
        column_title_side = "bottom", width = ncol(hcc_ps_orig_ident)*unit(4, "mm"), 
        height = nrow(hcc_ps_orig_ident)*unit(4, "mm"))
labAnch_hcc = FindTransferAnchors(reference = liver_and_cortex,query = hcc_sc_srat,dims=seq(30))
pred_hcc = TransferData(labAnch_hcc,refdata=liver_and_cortex@meta.data$annot,dims=seq(30))
hcc_sc_srat = AddMetaData(hcc_sc_srat,metadata=pred_hcc)

hcc_genes=read.table('/lustre/scratch117/casm/team274/gk14/Dediff/ho_hcc_sc/SRR14424777/output/Gene/filtered/features.tsv.gz', sep = "\t")

FeaturePlot(hcc_sc_srat,c('prediction.score.Adult_Hepatocytes', 'prediction.score.Adult_HHyP'))

FeaturePlot(hcc_sc_srat, c("EPCAM", "KRT19", "ALDH1A1", "ANPEP", "CD24", "CD44", "CD47", "THY1", "PROM1", "AFP"))
FeaturePlot(hcc_sc_srat, lit_genes)

DimPlot(hcc_sc_srat, group.by = "predicted.id")
DimPlot(hcc_sc_srat, group.by = "predicted.id", label = T)
DimPlot(hcc_sc_srat, label = T)
FeaturePlot(hcc_sc_srat, c("nFeature_RNA", "nCount_RNA", "percent.mt"))

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

saveRDS(hcc_sc_srat, '/lustre/scratch117/casm/team274/gk14/Dediff/hcc_sc.rds')
hcc_sc_srat=readRDS('/lustre/scratch117/casm/team274/gk14/Dediff/hcc_sc.rds')
DimPlot(hcc_sc_srat, group.by = "annot2")
colfunc <- colorRampPalette(c("#084594", "#C6DBEF"))
pdf("hcc_sc_umap.pdf", height = 3, width=3,useDingbats = F)
gg=DimPlot(hcc_sc_srat, group.by = "annot2", label=T,
           cols = c( colfunc(8),"#771122","#DDCC77"),
           pt.size = 0.1) +theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
                                 axis.ticks = element_blank(), legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1") +ggtitle("")
plot(gg)
dev.off()

DimPlot(hcc_sc_srat, group.by = "annot2", label=T,
        cols = c( colfunc(17),"#771122","#DDCC77"),
        pt.size = 0.1) +theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
                              axis.ticks = element_blank(), legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1") +ggtitle("")

DimPlot(hcc_sc_srat, group.by = "annot3", label=T,
        cols = c( colfunc(8),"#771122","#DDCC77"),
        pt.size = 0.1) +theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
                              axis.ticks = element_blank(), legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1") +ggtitle("")



inh_ps_bp=predictSimilarity(tr_inh_nomin, inhouse_hb_srat@assays$RNA@counts[inhouse_shared_genes,], logits = F)

hcc_ps_bp=predictSimilarity(train_for_hcc, hcc_sc_srat@assays$RNA@counts[hcc_inter,], logits = F)
canc_cells_hb=as.data.frame(inh_ps_bp[rownames(inhouse_hb_srat@meta.data)[which(inhouse_hb_srat@meta.data$annot=="Cancer cells")],])
canc_cells_hb$cell=rownames(canc_cells_hb)
canc_cells_hb$dataset="Hepatoblastoma cancer cells"
canc_cells_hb=canc_cells_hb[,c("cell","dataset", "Adult_Hepatocytes", "Fetal_Hepatocyte")]
canc_cells_hb$Adult_Hepatocytes=canc_cells_hb$Adult_Hepatocytes-0.5
canc_cells_hb$Fetal_Hepatocyte=canc_cells_hb$Fetal_Hepatocyte-0.5
canc_cells_hb$pos_or_neg_ahep="pos"
canc_cells_hb$pos_or_neg_ahep[which(canc_cells_hb$Adult_Hepatocytes<0)]="neg"
canc_cells_hb$pos_or_neg_fhep="pos"
canc_cells_hb$pos_or_neg_fhep[which(canc_cells_hb$Fetal_Hepatocyte<0)]="neg"

canc_cells_hcc=as.data.frame(hcc_ps_bp[rownames(hcc_sc_srat@meta.data)[which(hcc_sc_srat@meta.data$annot=="cancer")],])
canc_cells_hcc$cell=rownames(canc_cells_hcc)
canc_cells_hcc$dataset="Hepatocellular carcinoma cancer cells"
canc_cells_hcc=canc_cells_hcc[,c("cell","dataset", "Adult_Hepatocytes", "Fetal_Hepatocyte")]
canc_cells_hcc$Adult_Hepatocytes=canc_cells_hcc$Adult_Hepatocytes-0.5
canc_cells_hcc$Fetal_Hepatocyte=canc_cells_hcc$Fetal_Hepatocyte-0.5
canc_cells_hcc$pos_or_neg_ahep="pos"
canc_cells_hcc$pos_or_neg_ahep[which(canc_cells_hcc$Adult_Hepatocytes<0)]="neg"
canc_cells_hcc$pos_or_neg_fhep="pos"
canc_cells_hcc$pos_or_neg_fhep[which(canc_cells_hcc$Fetal_Hepatocyte<0)]="neg"




canc_cells_both=rbind(canc_cells_hcc, canc_cells_hb)
canc_cells_both_melted=melt(canc_cells_both, id.vars = c("cell", "dataset"))
canc_cells_both_melted$pos_or_neg="pos"
canc_cells_both_melted$pos_or_neg[which(canc_cells_both_melted$value<0)]="neg"
ggplot(canc_cells_both_melted,aes(y=reorder(cell, value), x=value)) + 
  geom_bar(stat = "identity") +facet_grid(vars(variable), vars(dataset), scales="free")+
  theme(text = element_text(size=12),axis.text.y = element_blank())

pdf('hb_fhep.pdf', height = 2, width = 2)
gg=ggplot(canc_cells_hb,aes(y=reorder(cell, Fetal_Hepatocyte), x=Fetal_Hepatocyte, color=pos_or_neg_fhep)) + 
  geom_bar(stat = "identity") +
  theme(text = element_text(size=12),axis.text.y = element_blank(), legend.position = "none" ) +xlim(c(-0.5,0.5))  + 
  scale_color_manual(values=c("#313695", "#a50026")) + ylab("HB cancer cells")
plot(gg)
dev.off()
pdf('hb_ahep.pdf', height = 2, width = 2)
gg=ggplot(canc_cells_hb,aes(y=reorder(cell, Adult_Hepatocytes), x=Adult_Hepatocytes, color=pos_or_neg_ahep)) + 
  geom_bar(stat = "identity") +
  theme(text = element_text(size=12),axis.text.y = element_blank(), legend.position = "none" ) +xlim(c(-0.5,0.5))  + 
  scale_color_manual(values=c("#313695", "#a50026")) + ylab("HB cancer cells")
plot(gg)
dev.off()
pdf('hcc_fhep.pdf', height = 2, width = 2)
gg=ggplot(canc_cells_hcc,aes(y=reorder(cell, Fetal_Hepatocyte), x=Fetal_Hepatocyte, color=pos_or_neg_fhep)) + 
  geom_bar(stat = "identity") +
  theme(text = element_text(size=12),axis.text.y = element_blank(), legend.position = "none" ) +xlim(c(-0.5,0.5))  + 
  scale_color_manual(values=c("#313695", "#a50026")) + ylab("HCC cancer cells")
plot(gg)
dev.off()
pdf('hcc_ahep.pdf', height = 2, width = 2)
gg=ggplot(canc_cells_hcc,aes(y=reorder(cell, Adult_Hepatocytes), x=Adult_Hepatocytes, color=pos_or_neg_ahep)) + 
  geom_bar(stat = "identity") +
  theme(text = element_text(size=12),axis.text.y = element_blank(), legend.position = "none" ) +xlim(c(-0.5,0.5))  + 
  scale_color_manual(values=c("#313695", "#a50026")) + ylab("HCC cancer cells")
plot(gg)
dev.off()

tr_inh1=trainModel(liver_and_cortex@assays$RNA@counts,
                  as.character(liver_and_cortex@meta.data$annot),
                  workers=NULL, minCells = 2000)
tr_inh2=trainModel(liver_and_cortex@assays$RNA@counts,
                   as.character(liver_and_cortex@meta.data$annot),
                   workers=NULL, minCells = 2000)
tr_inh3=trainModel(liver_and_cortex@assays$RNA@counts,
                   as.character(liver_and_cortex@meta.data$annot),
                   workers=NULL, minCells = 2000)


hcc_mtx=hcc_sc_srat@assays$RNA@counts
rownames(hcc_mtx)=hcc_genes$V1
adr_genes=ensembl_adrenal=read.table("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/SCP/Data/fAdrenal19wk/cellranger302_count_32644_WSSS_F_Adr8710632_GRCh38-1_2_0/filtered_feature_bc_matrix/features.tsv.gz", sep = "\t", header = F)
adr_genes_l=adr_genes
rownames(adr_genes_l)=make.unique(adr_genes_l$V2)
adr_genes_l=adr_genes_l[rownames(liver_gastr_ee),]
rownames(adr_genes_l)=adr_genes_l$V1
inter_hcc_landc=intersect(rownames(hcc_mtx), rownames(adr_genes_l))
hcc_mtx=hcc_mtx[inter_hcc_landc,]
adr_genes_l_hcc=adr_genes_l[inter_hcc_landc,]
identical(rownames(hcc_mtx), rownames(adr_genes_l_hcc))

rownames(hcc_mtx)=adr_genes_l_hcc$V2

ps_n2000=predictSimilarity(tr_inh1, hcc_mtx, minGeneMatch = 0.7)
similarityHeatmap(ps_n2000, row_split=hcc_sc_srat@meta.data$annot3, use_raster=F, row_title_rot=0)


pr_inh=predictSimilarity(tr_inh_ee_gast, hcc_mtx, minGeneMatch = 0.9)

similarityHeatmap(pr_inh, row_split=hcc_sc_srat@meta.data$annot3, use_raster=F, row_title_rot=0)

sort(table(liver_gastr_ee@meta.data$annot))

liver_hep_only=subset(liver_gastr_ee, annot%in%c("Fetal_Hepatocyte", "Adult_Hepatocytes"))
not_in=unique(liver_gastr_ee@meta.data$annot[!(liver_gastr_ee@meta.data$annot%in%c("Fetal_Hepatocyte", "Adult_Hepatocytes"))])
liver_no_hep=subset(liver_gastr_ee, annot%in%not_in)
Idents(liver_no_hep)=liver_no_hep@meta.data$annot
liver_no_hep_downs=subset(liver_no_hep, downsample=100)
downs_and_hep=merge(liver_no_hep_downs, liver_hep_only)

tr1=trainModel(downs_and_hep@assays$RNA@counts, 
               as.character(downs_and_hep@meta.data$annot),
               workers = NULL, minCells = Inf)

tr2=trainModel(downs_and_hep@assays$RNA@counts, 
               as.character(downs_and_hep@meta.data$annot),
               workers = NULL, minCells = Inf)

tr3=trainModel(downs_and_hep@assays$RNA@counts, 
               as.character(downs_and_hep@meta.data$annot),
               workers = NULL, minCells = Inf)
pr_inh=predictSimilarity(tr1, hcc_mtx, hcc_sc_srat@meta.data$annot3, minGeneMatch = 0.9,logits = F)
pr_in2=predictSimilarity(tr2, hcc_mtx, hcc_sc_srat@meta.data$annot3, minGeneMatch = 0.9, logits = F)
pr_in3=predictSimilarity(tr3, hcc_mtx, hcc_sc_srat@meta.data$annot3, minGeneMatch = 0.9, logits = F)

sg=intersect(rownames(inhouse_hb_srat), rownames(downs_and_hep))
pr_hb=predictSimilarity(tr1, inhouse_hb_srat@assays$RNA@counts[sg,], inhouse_hb_srat@meta.data$annot, logits = F)

similarityHeatmap(pr_inh)
similarityHeatmap(pr_in2)
similarityHeatmap(pr_in3)
similarityHeatmap(pr_hb)

col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
Heatmap(pr_inh, cluster_columns=F, cluster_rows = F, col=col_fun, 
        column_gap = unit(3, "mm"), row_names_side = "left",
        column_title_side = "bottom", width = ncol(pr_inh)*unit(4, "mm"), 
        height = nrow(pr_inh)*unit(4, "mm"))

Heatmap(pr_hb, cluster_columns=F, cluster_rows = F, col=col_fun, 
        column_gap = unit(3, "mm"), row_names_side = "left",
        column_title_side = "bottom", width = ncol(pr_hb)*unit(4, "mm"), 
        height = nrow(pr_hb)*unit(4, "mm"))

train_liver=trainModel(liver_gastr_ee@assays$RNA@counts, liver_gastr_ee@meta.data$annot, workers = NULL)

pr_hcc=predictSimilarity(train_liver, hcc_sc_srat@assays$RNA@counts, hcc_sc_srat@meta.data$annot3, minGeneMatch = 0.8, logits = F)
similarityHeatmap(pr_hcc)
pr_inh=predictSimilarity(train_liver, inhouse_hb_srat@assays$RNA@counts, inhouse_hb_srat@meta.data$annot, minGeneMatch = 0.4, logits=F)

similarityHeatmap(pr_inh)

Heatmap(pr_hcc, cluster_columns=F, cluster_rows = F, column_split=colsplit_df$wide_annot, col=col_fun, 
        column_order=levels(colsplit_df$celltype), column_gap = unit(3, "mm"), row_names_side = "left",
        column_title_side = "bottom", width = ncol(pr_hcc)*unit(4, "mm"), 
        height = nrow(pr_hcc)*unit(4, "mm"))

Heatmap(pr_inh, cluster_columns=F, cluster_rows = F, column_split=colsplit_df$wide_annot, col=col_fun, 
        column_order=levels(colsplit_df$celltype), column_gap = unit(3, "mm"), row_names_side = "left",
        column_title_side = "bottom", width = ncol(pr_inh)*unit(4, "mm"), 
        height = nrow(pr_inh)*unit(4, "mm"))

colsplit_df=data.frame(celltype=colnames(pr_inh))
colsplit_df$wide_annot=colsplit_df$celltype
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%c("Fetal_Mono-Mac","Fetal_B-cells","Fetal_DCs",
                                                         "Fetal_Monocyte",
                                                         "Fetal_Early lymphoid_T lymphocyte",
                                                         "Fetal_VCAM1+ EI macrophage", "Fetal_ILC precursor",
                                                         "Fetal_Neutrophil-myeloid progenitor",
                                                         "Fetal_Monocyte precursor", 
                                                         "Fetal_NK",
                                                         "Fetal_Mast cell", 
                                                         "Fetal_Kupffer Cell"))]="Fetal_leukocytes"
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%c("Adult_Non-inflammatory_macs", 
                                                         "Adult_T-cells", 
                                                         "Adult_Inflammatory_macs",
                                                         "Adult_NK-like_cells",
                                                         "Adult_Mature_Bcells"))]="Adult_leukocytes"
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%c("Adult_Hepatocytes", 
                                                         "Fetal_Hepatocyte", 
                                                         "Adult_HHyP",
                                                         "Fetal_HHyP",
                                                         "Adult_BECs"))]="EPCAM+"

colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%c("Fetal_Erythroid", 
                                                         "Fetal_MEMP", 
                                                         "Fetal_Megakaryocyte",
                                                         "Fetal_MEMP",
                                                         "Fetal_HSC_MPP",
                                                         "Adult_Plasma_cells",
                                                         "Adult_Erythroid_cells"))]="Other blood cells"
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%c("Cortex"))]="Negative control"
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%c("Hypoblast", "Epiblast"))]="Early embryo"
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%c("Mesoderm", "Endoderm", "Non-Neural Ectoderm"))]="Gastrulation"
colsplit_df$wide_annot[which(colsplit_df$wide_annot%in%c("Fetal_Endothelial cell",
                                                         "Adult_LSECs",
                                                         "Adult_Portal_endothelium"))]="Endothelium"

colsplit_df$wide_annot=factor(colsplit_df$wide_annot, levels=c("Early embryo", "Gastrulation","EPCAM+", "Adult_Hepatic_stellate_cells", "Endothelium",
                                                               "Other blood cells","Adult_leukocytes", "Fetal_leukocytes", 
                                                               "Fetal_Fibroblast",  "Negative control" ))
colsplit_df$celltype=factor(colsplit_df$celltype, levels = c("Hypoblast", "Epiblast","Mesoderm", "Endoderm", "Non-Neural Ectoderm","Fetal_Hepatocyte",
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
                                                             "Cortex"))


train_liver1=trainModel(liver_gastr_ee@assays$RNA@counts, liver_gastr_ee@meta.data$annot, workers = NULL)
train_liver2=trainModel(liver_gastr_ee@assays$RNA@counts, liver_gastr_ee@meta.data$annot, workers = NULL)
train_liver3=trainModel(liver_gastr_ee@assays$RNA@counts, liver_gastr_ee@meta.data$annot, workers = NULL)

pr_hcc1=predictSimilarity(train_liver1, hcc_sc_srat@assays$RNA@counts, hcc_sc_srat@meta.data$annot3, minGeneMatch = 0.8, logits = F)
pr_hcc2=predictSimilarity(train_liver2, hcc_sc_srat@assays$RNA@counts, hcc_sc_srat@meta.data$annot3, minGeneMatch = 0.8, logits = F)
pr_hcc3=predictSimilarity(train_liver3, hcc_sc_srat@assays$RNA@counts, hcc_sc_srat@meta.data$annot3, minGeneMatch = 0.8, logits = F)
pr_inh1=predictSimilarity(train_liver1, inhouse_hb_srat@assays$RNA@counts, inhouse_hb_srat@meta.data$annot, minGeneMatch = 0.4, logits=F)
pr_inh2=predictSimilarity(train_liver2, inhouse_hb_srat@assays$RNA@counts, inhouse_hb_srat@meta.data$annot, minGeneMatch = 0.4, logits=F)
pr_inh3=predictSimilarity(train_liver3, inhouse_hb_srat@assays$RNA@counts, inhouse_hb_srat@meta.data$annot, minGeneMatch = 0.4, logits=F)

Heatmap(pr_hcc1[,as.character(colsplit_df$celltype)], cluster_columns=F, cluster_rows = F, column_split=colsplit_df$wide_annot, col=col_fun, 
        column_order=levels(colsplit_df$celltype), column_gap = unit(3, "mm"), row_names_side = "left",
        column_title_side = "bottom", width = ncol(pr_hcc1)*unit(4, "mm"), 
        height = nrow(pr_hcc1)*unit(4, "mm"))

Heatmap(pr_inh1[,as.character(colsplit_df$celltype)], cluster_columns=F, cluster_rows = F, column_split=colsplit_df$wide_annot, col=col_fun, 
        column_order=levels(colsplit_df$celltype), column_gap = unit(3, "mm"), row_names_side = "left",
        column_title_side = "bottom", width = ncol(pr_inh)*unit(4, "mm"), 
        height = nrow(pr_inh)*unit(4, "mm"))

Heatmap(pr_hcc2[,as.character(colsplit_df$celltype)], cluster_columns=F, cluster_rows = F, column_split=colsplit_df$wide_annot, col=col_fun, 
        column_order=levels(colsplit_df$celltype), column_gap = unit(3, "mm"), row_names_side = "left",
        column_title_side = "bottom", width = ncol(pr_hcc1)*unit(4, "mm"), 
        height = nrow(pr_hcc1)*unit(4, "mm"))

Heatmap(pr_inh2[,as.character(colsplit_df$celltype)], cluster_columns=F, cluster_rows = F, column_split=colsplit_df$wide_annot, col=col_fun, 
        column_order=levels(colsplit_df$celltype), column_gap = unit(3, "mm"), row_names_side = "left",
        column_title_side = "bottom", width = ncol(pr_inh)*unit(4, "mm"), 
        height = nrow(pr_inh)*unit(4, "mm"))

Heatmap(pr_hcc3[,as.character(colsplit_df$celltype)], cluster_columns=F, cluster_rows = F, column_split=colsplit_df$wide_annot, col=col_fun, 
        column_order=levels(colsplit_df$celltype), column_gap = unit(3, "mm"), row_names_side = "left",
        column_title_side = "bottom", width = ncol(pr_hcc1)*unit(4, "mm"), 
        height = nrow(pr_hcc1)*unit(4, "mm"))

Heatmap(pr_inh3[,as.character(colsplit_df$celltype)], cluster_columns=F, cluster_rows = F, column_split=colsplit_df$wide_annot, col=col_fun, 
        column_order=levels(colsplit_df$celltype), column_gap = unit(3, "mm"), row_names_side = "left",
        column_title_side = "bottom", width = ncol(pr_inh)*unit(4, "mm"), 
        height = nrow(pr_inh)*unit(4, "mm"))

similarityHeatmap(pr_inh3)

