###############################
# tabula sapiens analysis #####
##############################

message('now install remotes')
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "http://cran.us.r-project.org")
}
remotes::install_github("mojaveazure/seurat-disk")
message('installed seurat-disk')
install.packages('remotes')
remotes::install_github(repo = 'satijalab/seurat', ref = 'develop')
message('installed seurat')

library(SeuratDisk)
library(Seurat)

Convert("~/TabulaSapiens.h5ad", "~/TabulaSapiens.h5seurat")
message("Tabula converted")
# This creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory

# This .d5seurat object can then be read in manually
seuratObject <- LoadH5Seurat("~/TabulaSapiens.h5seurat", assays='RNA')
write.table(seuratObject@meta.data,"TS_metadata.tsv",row.names=TRUE,col.names=TRUE,sep='\t',quote=FALSE )
message("Tabula loaded & metadata written")
#seuratObject@meta.data$free_annotation
cellSigs = data.frame(do.call(cbind,lapply(split(rownames(seuratObject@meta.data),seuratObject@meta.data$free_annotation),function(e) rowSums(seuratObject@assays$RNA@counts[,e,drop=FALSE]))))
rm(seuratObject) # not necessary but to save memory
gns <- read.table("~/key_matching_gene_names_with_ensids.txt",header=T)
message("gns read")
#intersect
cellSigs<-cellSigs[which(rownames(cellSigs) %in% gns$SYMBOL),]
idx<-match(rownames(cellSigs),gns$SYMBOL)
rownames(cellSigs)<-gns$GENEID[idx]
write.table(cbind(geneName=rownames(cellSigs),cellSigs),'/nfs/research/icortes/mkalyva/cellSigs_summary_TS_sapiens_without_norm.tsv',row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
# normalise atlas reference to 1
cellSigs = as.data.frame(t(t(cellSigs)/colSums(cellSigs)))
message("reference ready")
write.table(cbind(geneName=rownames(cellSigs),cellSigs),'/nfs/research/icortes/mkalyva/cellSigs_summary_TS_sapiens.tsv',row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
message("done")

