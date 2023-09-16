library(argparse)
library(Seurat)
# make reference for cibersort

args = commandArgs(trailingOnly=TRUE)

data<-read.csv(args[1],header=T,stringsAsFactors = F, sep=",") #  has to be an "rmbatchdge.txt" as downloaded from the paper's figshare 
message("data read")
ann<-read.csv(args[2],header=T,stringsAsFactors = F, sep=",") # has to be an "rmbatchAnno.csv" as downloaded from the paper's figshare 
message("ann read")
name<-ann$Batch[1]
rownames(ann)<-ann$Cell_id
seuratObject<-CreateSeuratObject(counts = data,meta.data = ann)
message("seurat done read")
cellSigs = data.frame(do.call(cbind,lapply(split(rownames(seuratObject@meta.data),seuratObject@meta.data$Celltype),function(e) rowSums(seuratObject@assays$RNA@counts[,e,drop=FALSE]))))
message("cellSigs done")
gns <- read.table("~/key_matching_gene_names_with_ensids.txt",header=T, sep="\t")
message("gns read")
#intersect
cellSigs<-cellSigs[which(rownames(cellSigs) %in% gns$SYMBOL),]
idx<-match(rownames(cellSigs),gns$SYMBOL)
rownames(cellSigs)<-gns$GENEID[idx]
write.table(cbind(geneName=rownames(cellSigs),cellSigs),paste0("~/Han_Nature_2020/cellSigs/cellSigs_before_norm",name,".txt"),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
# normalise to 1
cellSigs = as.data.frame(t(t(cellSigs)/colSums(cellSigs)))
message("reference ready")
write.table(cbind(geneName=rownames(cellSigs),cellSigs),paste0("~Han_Nature_2020/cellSigs/cellSigs_",name,".txt"),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
message("done")