#' Calculates the stemness score on a range of samples to compare to the deconvolution results
setwd('~/trueHome/Projects/stemness/')

#############
# Libraries #
#############

library(biomaRt)
library(tximport)
library(EnsDb.Hsapiens.v86)

##########
# Params #
##########

stemness = read.table('Data/stemnessGeneWeights.tsv',sep='\t',header=FALSE)
stemness = setNames(stemness[,2],stemness[,1])
bigBulkMani = readRDS('Data/bigBulkManifest.RDS')
stJudesSrc = 'Data/StJudesTxTPMs.tsv'
stJudesMeta = 'Data/StJudesBasicMetadata.tsv'

#############
# Functions #
#############

#' Load bulk data
#'
#' Load bulk data from a list of samples stored in the required format.
#'
#' @param tgts Files to load.
#' @return A list containing counts and lengths.
loadBulk = function(tgts){
  #Load them
  dat = lapply(tgts,read.table,sep='\t',header=TRUE)
  gns = Reduce(intersect,lapply(dat,function(e) e$geneName))
  cnts = do.call(cbind,lapply(dat,function(e) e[match(gns,e$geneName),-c(1,2),drop=FALSE])) 
  rownames(cnts) = gns
  lens = do.call(cbind,lapply(dat,function(e) {
                                tmp = e[match(gns,e$geneName),rep(match('geneLengths',colnames(e)),ncol(e)-2),drop=FALSE] 
                                colnames(tmp) = colnames(e)[-c(1,2)]
                                tmp}))
  rownames(lens) = gns
  return(list(cnts=cnts,lens=lens))
}

#' Calculates stemness from bulk transcriptomes
#'
#' Uses the supplied weights to calculate stemness.  Note that given the definition of stemness from that paper, this is a relative, not absolute, measure.  As such you will get different results depending on what samples you run together.
#'
#' @param tpm A matrix with genes as rows and samples as columns giving gene expression in  TPM format.
#' @param stemnessGeneWeights The gene weights to use to calculate stemness.
#' @param convertGeneNames If true, assumes gene names are given in ENSEMBL format and converts to the assumed gene symbol naming of the stemness weights.
#' @return A vector of stemness scores per sample.
calcStemness = function(tpm,stemnessGeneWeights=stemness,convertGeneNames=TRUE){
  genes = rownames(tpm)
  nameConv = getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol",'description'),values=genes,mart= mart)
  #Munge stemness score
  m = match(names(stemness),nameConv$hgnc_symbol)
  names(stemness) = nameConv$ensembl_gene_id[m]
  stemness = stemness[!is.na(names(stemness))]
  #Calculate the stemness things on the bulk data
  tmp = tpm[names(stemness),]
  tmp = scale(tmp)
  x = cor(tmp,stemness,method='sp')
  x = (x-min(x))/(max(x)-min(x))
  return(x[,1])
}



######################
# Calculate Stemness #
######################

#Load adult data and convert to tpms
mDat = bigBulkMani[bigBulkMani$biopsy=='Primary Tumor' & bigBulkMani$dataset=='TCGA',]
bulkDat = loadBulk(mDat$fullPath)
bulkDat = bulkDat$cnts/bulkDat$lens
bulkDat = t(t(bulkDat)/colSums(bulkDat))*1e9
#Load StJudes and convert to gene level.  We've lost the intermediate files need for tximport to do this for us, so do it manually
sjDat = read.delim(stJudesSrc,sep='\t',header=TRUE)
tx2gene = transcripts(EnsDb.Hsapiens.v86,return.type='DataFrame')
rownames(sjDat) = gsub('\\.[0-9]*$','',sjDat[,1])
sjDat = sjDat[,-1]
tx2gene = tx2gene[,c('tx_id','gene_id')]
colnames(tx2gene) = c('TXNAME','GENEID')
txToKeep = rownames(sjDat) %in% tx2gene$TXNAME
sjDat = sjDat[txToKeep,]
m = match(rownames(sjDat),tx2gene$TXNAME)
sjDat = rowsum(sjDat,tx2gene$GENEID[m])
#Merge and log transform
gns = intersect(rownames(bulkDat),rownames(sjDat))
tpm = cbind(bulkDat[gns,],sjDat[gns,])
tpm = log10(1+tpm)
#Get stemness
stemScore = calcStemness(tpm)
#Get stjudes metadata
sjMeta = read.delim(stJudesMeta,sep='\t',header=TRUE)
sjMeta = sjMeta[sjMeta$sequencing_type=='RNA-SEQ' & grepl('\\.bam$',sjMeta$file_path),]
m = match(gsub('\\..*','',colnames(sjDat)),sjMeta$sample_name)
sjMeta = sjMeta[m,]
#Determine groupings
groups = c(mDat$disease,sjMeta$attr_diagnosis)
tmp = split(stemScore,groups)
tmp = tmp[order(sapply(tmp,median))]
boxplot(tmp,las=2,col=ifelse(names(tmp) %in% mDat$disease,'white','darkred'),horizontal=TRUE)
dd = data.frame(type = ifelse(names(tmp) %in% mDat$disease,'TCGA','StJudes'),disease=names(tmp),score = sapply(tmp,median))
rownames(dd) = NULL










