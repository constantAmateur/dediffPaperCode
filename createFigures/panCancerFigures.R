#' Make figures analysing pan-cancer data
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

#############
# Functions #
#############

#' Calculates stemness from bulk transcriptomes
#'
#' Uses the supplied weights to calculate stemness.  Note that given the definition of stemness from that paper, this is a relative, not absolute, measure.  As such you will get different results depending on what samples you run together.
#'
#' @param tpm A matrix with genes as rows and samples as columns giving gene expression in  TPM format.
#' @param stemnessGeneWeights The gene weights to use to calculate stemness.
#' @param convertGeneNames If true, assumes gene names are given in ENSEMBL format and converts to the assumed gene symbol naming of the stemness weights.
#' @return A vector of stemness scores per sample.
calcStemness = function(tpm,stemnessGeneWeights=stemness,convertGeneNames=TRUE){
  if(convertGeneNames){
    genes = rownames(tpm)
    nameConv = getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol",'description'),values=genes,mart= mart)
    #Munge stemness score
    m = match(names(stemness),nameConv$hgnc_symbol)
    names(stemness) = nameConv$ensembl_gene_id[m]
    stemness = stemness[!is.na(names(stemness))]
    #Calculate the stemness things on the bulk data
    tmp = tpm[names(stemness),]
  }else{
    gns = intersect(rownames(tpm),names(stemness))
    tmp = tpm[gns,]
    stemness = stemness[gns]
  }
  tmp = scale(tmp)
  x = cor(tmp,stemness,method='sp')
  x = (x-min(x))/(max(x)-min(x))
  return(x[,1])
}



##########
# Params #
##########

nParallel=48
plotDir = 'Results/panFigures/'
srcCS = 'Results/cellSignalOutput/Marias/'

#############
# pR2 vs EE #
#############

#Load the TCGA data
dat = list.files(file.path(srcCS,'TCGA','results_TCGA'))
dat = grep('EE_',dat,invert=TRUE,value=TRUE)
#Get the cancer types
cTypes = unique(gsub('_(Fetal|TS|ee)_.*','',dat))
outs = list()
toPlot = c('pR2','ee')
for(cType in cTypes){
  message(cType)
  #Load early embryo
  ee = csa$normaliseExposures(file.path(srcCS,'TCGA','results_TCGA',paste0(cType,'_ee_fitExposures.tsv')))
  ee = data.frame(sample = colnames(ee$exposures),
                  ee = colSums(ee$exposures[c('Epiblast','Hypoblast'),]),
                  int = as.numeric(ee$exposures["Intercept",]),
                  pR2 = as.numeric(ee$gof['pR2',]))
  #Add gastrulation
  gast = csa$normaliseExposures(file.path(srcCS,'TCGA','results_TCGA',paste0(cType,'_ee_gast_fitExposures.tsv')))
  gast = data.frame(sample = ee$sample,
                    #ee = colSums(gast$exposures[c("Mesoderm","Epiblast", "Hypoblast"),ee$sample]),
                    gast = colSums(gast$exposures[c("Endoderm", "Mesoderm", "Non-Neural Ectoderm", "Epiblast", "Hypoblast"),ee$sample]),
                    ee = colSums(gast$exposures[c('Epiblast','Hypoblast'),ee$sample]),
                    int = as.numeric(gast$exposures["Intercept",ee$sample]),
                    pR2 = as.numeric(gast$gof['pR2',ee$sample]))
  #Full
  adult = csa$normaliseExposures(file.path(srcCS,'TCGA','results_TCGA',paste0(cType,'_TS_fitExposures.tsv')))
  adult = data.frame(sample = ee$sample,
                     #ee = colSums(adult$exposures[c("Mesoderm","Epiblast", "Hypoblast"),ee$sample]),
                     gast = colSums(adult$exposures[c("Endoderm", "Mesoderm","Non.NeuralEctoderm", "Epiblast", "Hypoblast"),ee$sample]),
                     ee = colSums(adult$exposures[c('Epiblast','Hypoblast'),ee$sample]),
                     int = as.numeric(adult$exposures["Intercept",ee$sample]),
                     pR2 = as.numeric(adult$gof['pR2',ee$sample]))
  outs[[cType]] = list(ee=ee,gast=gast,adult=adult)
}
#Get numbers for manuscript
t1 = do.call(rbind,lapply(outs[!names(outs) %in% 'HSC_blast'],function(e) e$ee))
#mean(t1$pR2)
#mean(t1$int)
#Custom is my version where the adult shows the gastrulation signal not the EE and we only show endpoints for the blast.  standard is standard
plotType='custom'
#Plot the arrows?
plotArrows=FALSE
#Some other params
specType='HSC_blast'
cols = list(ee='#4DBBD5',
            gast='#8F3931',
            adult='#CC79A7',
            arrow='#000000',
            arrowBlast='#8F3931')
#855C75,#D9AF6B,#AF6458,#736F4C,#526A83,#625377,#68855C,#9C9C5E,#A06177,#8C785D,#467378,#7C7C7C
cols$ee = '#D9AF6B'
cols$gast = '#AF6458'
cols$adult = '#68855C'
cols$fetal = '#855C75'
alphas = list(arrow=as.hexmode(255*1),
              pnt=as.hexmode(255*0.2),
              pntBlast=as.hexmode(255*1))
lwds = list(arrow=0.5,
            arrowBlast=0.3)
cexs = list(pnt = 0.1,
            pntBlast=0.3)
arrLen=0.08
pdf(file.path(plotDir,'embryoSignal.pdf'),width=4,height=3)
par(mar=c(5,4,1,1))
if(plotType=='custom'){
  plot(NA,
       type='n',
       xlab='pR2',
       ylab='EE or Gastrulation',
       xlim=c(0,1),
       ylim=c(0,1))
  for(cType in cTypes){
    ee = outs[[cType]]$ee
    gast = outs[[cType]]$gast
    adult = outs[[cType]]$adult
    if(plotArrows){
      if(cType==specType){
        arrows(x0 = (ee[,'pR2']),
               y0 = (ee[,'ee']),
               x1 = (adult[,'pR2']),
               y1 = (adult[,'gast']),
               lwd=lwds$arrowBlast,
               col=paste0(cols$arrowBlast,alphas$arrow),
               length=arrLen)
        #arrows(x0 = (gast[,'pR2']),
        #       y0 = (gast[,'ee']),
        #       x1 = (adult[,'pR2']),
        #       y1 = (adult[,'gast']),
        #       lwd=lwds$arrowBlast,
        #       col=paste0(cols$arrowBlast,alphas$arrow),
        #       length=arrLen)
      }else{
        arrows(x0 = mean(ee[,'pR2']),
               y0 = mean(ee[,'ee']),
               x1 = mean(gast[,'pR2']),
               y1 = mean(gast[,'ee']),
               lwd=lwds$arrow,
               col=paste0(cols$arrow,alphas$arrow),
               length=arrLen)
        arrows(x0 = mean(gast[,'pR2']),
               y0 = mean(gast[,'ee']),
               x1 = mean(adult[,'pR2']),
               y1 = mean(adult[,'gast']),
               lwd=lwds$arrow,
               col=paste0(cols$arrow,alphas$arrow),
               length=arrLen)
      }
      #arrows(x0 = ee[,toPlot[1]],
      #       y0 = ee[,toPlot[2]],
      #       x1 = gast[,toPlot[1]],
      #       y1 = gast[,toPlot[2]],
      #       lwd=arrLwd,
      #       col=arrCol,
      #       length=arrLen)
      #arrows(x0 = gast[,toPlot[1]],
      #       y0 = gast[,toPlot[2]],
      #       x1 = adult[,toPlot[1]],
      #       y1 = adult[,toPlot[2]],
      #       lwd=arrLwd,
      #       col=arrCol,
      #       length=arrLen)
    }
  }
  #Plot points
  #for(cType in cTypes){
  #  points(outs[[cType]]$gast$pR2-outs[[cType]]$ee$pR2,(outs[[cType]]$gast$ee-outs[[cType]]$ee$ee)/outs[[cType]]$ee$ee,pch=19,cex=0.1,col='#1b9e7744')
  #  points(outs[[cType]]$adult$pR2-outs[[cType]]$gast$pR2,(outs[[cType]]$adult$gast-outs[[cType]]$gast$gast)/outs[[cType]]$gast$gast,pch=19,cex=0.1,col='#d95f0244')
  #}
  #cType='HSC_blast'
  #points(outs[[cType]]$gast$pR2-outs[[cType]]$ee$pR2,(outs[[cType]]$gast$ee-outs[[cType]]$ee$ee)/outs[[cType]]$ee$ee,pch=19,cex=0.5,col='red')
  #points(outs[[cType]]$adult$pR2-outs[[cType]]$gast$pR2,(outs[[cType]]$adult$gast-outs[[cType]]$gast$gast)/outs[[cType]]$gast$gast,pch=19,cex=0.5,col='green')
  specType='HSC_blast'
  for(cType in cTypes){
    if(cType==specType)
      next
    points(outs[[cType]]$ee[,c('pR2','ee')],pch=19,cex=cexs$pnt,col=paste0(cols$ee,alphas$pnt))
    points(outs[[cType]]$gast[,c('pR2','ee')],pch=19,cex=cexs$pnt,col=paste0(cols$gast,alphas$pnt))
    points(outs[[cType]]$adult[,c('pR2','gast')],pch=19,cex=cexs$pnt,col=paste0(cols$adult,alphas$pnt))
  }
  cType=specType
  points(outs[[cType]]$ee[,c('pR2','ee')],pch=19,cex=cexs$pntBlast,col=paste0(cols$ee,alphas$pntBlast))
  #points(outs[[cType]]$gast[,c('pR2','ee')],pch=19,cex=1,col='#d95f02FF')
  points(outs[[cType]]$adult[,c('pR2','gast')],pch=19,cex=cexs$pntBlast,col=paste0(cols$adult,alphas$pntBlast))
}else{
  plot(NA,
       type='n',
       xlab='pR2',
       ylab='EE',
       xlim=c(0,1),
       ylim=c(0,1))
  for(cType in cTypes){
    ee = outs[[cType]]$ee
    gast = outs[[cType]]$gast
    adult = outs[[cType]]$adult
    if(plotArrows){
      if(cType==specType){
        arrows(x0 = (ee[,'pR2']),
               y0 = (ee[,'ee']),
               x1 = (gast[,'pR2']),
               y1 = (gast[,'ee']),
               lwd=lwds$arrowBlast,
               col=paste0(cols$arrow,alphas$arrow),
               length=arrLen)
        arrows(x0 = (gast[,'pR2']),
               y0 = (gast[,'ee']),
               x1 = (adult[,'pR2']),
               y1 = (adult[,'ee']),
               lwd=lwds$arrowBlast,
               col=paste0(cols$arrow,alphas$arrow),
               length=arrLen)
      }else{
        arrows(x0 = mean(ee[,'pR2']),
               y0 = mean(ee[,'ee']),
               x1 = mean(gast[,'pR2']),
               y1 = mean(gast[,'ee']),
               lwd=lwds$arrow,
               col=paste0(cols$arrow,alphas$arrow),
               length=arrLen)
        arrows(x0 = mean(gast[,'pR2']),
               y0 = mean(gast[,'ee']),
               x1 = mean(adult[,'pR2']),
               y1 = mean(adult[,'ee']),
               lwd=lwds$arrow,
               col=paste0(cols$arrow,alphas$arrow),
               length=arrLen)
      }
      #arrows(x0 = ee[,toPlot[1]],
      #       y0 = ee[,toPlot[2]],
      #       x1 = gast[,toPlot[1]],
      #       y1 = gast[,toPlot[2]],
      #       lwd=arrLwd,
      #       col=arrCol,
      #       length=arrLen)
      #arrows(x0 = gast[,toPlot[1]],
      #       y0 = gast[,toPlot[2]],
      #       x1 = adult[,toPlot[1]],
      #       y1 = adult[,toPlot[2]],
      #       lwd=arrLwd,
      #       col=arrCol,
      #       length=arrLen)
    }
  }
  #Plot points
  #for(cType in cTypes){
  #  points(outs[[cType]]$gast$pR2-outs[[cType]]$ee$pR2,(outs[[cType]]$gast$ee-outs[[cType]]$ee$ee)/outs[[cType]]$ee$ee,pch=19,cex=0.1,col='#1b9e7744')
  #  points(outs[[cType]]$adult$pR2-outs[[cType]]$gast$pR2,(outs[[cType]]$adult$gast-outs[[cType]]$gast$gast)/outs[[cType]]$gast$gast,pch=19,cex=0.1,col='#d95f0244')
  #}
  #cType='HSC_blast'
  #points(outs[[cType]]$gast$pR2-outs[[cType]]$ee$pR2,(outs[[cType]]$gast$ee-outs[[cType]]$ee$ee)/outs[[cType]]$ee$ee,pch=19,cex=0.5,col='red')
  #points(outs[[cType]]$adult$pR2-outs[[cType]]$gast$pR2,(outs[[cType]]$adult$gast-outs[[cType]]$gast$gast)/outs[[cType]]$gast$gast,pch=19,cex=0.5,col='green')
  specType='HSC_blast'
  for(cType in cTypes){
    if(cType==specType)
      next
    points(outs[[cType]]$ee[,c('pR2','ee')],pch=19,cex=cexs$pnt,col=paste0(cols$ee,alphas$pnt))
    points(outs[[cType]]$gast[,c('pR2','ee')],pch=19,cex=cexs$pnt,col=paste0(cols$gast,alphas$pnt))
    points(outs[[cType]]$adult[,c('pR2','ee')],pch=19,cex=cexs$pnt,col=paste0(cols$adult,alphas$pnt))
  }
  cType=specType
  points(outs[[cType]]$ee[,c('pR2','ee')],pch=19,cex=cexs$pntBlast,col=paste0(cols$ee,alphas$pntBlast))
  points(outs[[cType]]$gast[,c('pR2','ee')],pch=19,cex=cexs$pntBlast,col=paste0(cols$gast,alphas$pntBlast))
  points(outs[[cType]]$adult[,c('pR2','ee')],pch=19,cex=cexs$pntBlast,col=paste0(cols$adult,alphas$pntBlast))
}
dev.off()

################################
# Childhood vs adult fetalness #
################################

#################
# TARGET cleanup
target = csa$normaliseExposures(file.path(srcCS,'PEDIATRIC','pediatric_results','PED_TS_Fetal_new_fitExposures.tsv'))
#Drop those that are not primary tumours
x=colnames(target$exposures)
tmp = strsplit(x,'\\.')
isBad = rep(FALSE,length(x))
#Drop those where there's two samples in one
isBad[lengths(tmp)>8]=TRUE
#Drop those that aren't TARGET
isBad[!grepl("^TARGET",x)]=TRUE
#Drop normal, AML,
isBad[which(sapply(tmp,function(e) e[3]) %in% c('00','20'))]=TRUE
#Drop bad tissue codes: recurrence, blood cancer, Mets, normal, cell line,
isBad[which(sapply(tmp,function(e) e[5])!='01A')]=TRUE
#Now we have to drop the many(!?) duplicates
isBad[duplicated(sapply(tmp,function(e) paste0(e[4],'.',e[5])))]=TRUE
#Finish, keep the good ones
goodSamps = x[!isBad]
#Make some basic meta-data
mDat = data.frame(sampleID = goodSamps,
                  disease = c(NBL='Neuroblastoma',RT='Rhabdoid Tumor',WT='Wilms Tumor')[gsub('TARGET\\.([A-Z]+)_.*','\\1',goodSamps)])
rownames(mDat) = mDat$sampleID
#Finalise
target$exposures = target$exposures[,goodSamps]
target$gof = target$gof[,goodSamps]
target$raw = target$raw[,goodSamps]
target$mDat = mDat
##################
# StJudes cleanup
judes = csa$normaliseExposures(file.path(srcCS,'PEDIATRIC','pediatric_results','Judes_TS_Fetal_new_fitExposures.tsv'))
#Get meta-data
mDat = read.delim('Data/StJudesBasicMetadata.tsv',sep='\t',header=TRUE)
#Drop the columns that are not related to RNA-seq
mDat = mDat[mDat$sequencing_type=='RNA-SEQ',]
#Drop the bam index duplicate rows
mDat = mDat[gsub('.*\\.','',mDat$file_path)=='bam',]
#There are 23 duplicate entries, seem to just be the same patient/data in two datasets.  Drop the smaller dataset entry
mDat = mDat[!(mDat$sample_name %in% mDat$sample_name[duplicated(mDat$sample_name)] & mDat$sj_datasets=='Clinical Pilot'),]
#We now have 1 row per sample_name, keep just our samples
x=colnames(judes$exposures)
mDat = mDat[match(x,mDat$sample_name),]
#Drop those that are not diagnostic samples
mDat = mDat[mDat$sample_type=='Diagnosis',]
#Get codes.  Save to prevent instability by fetching each time
tgtFile = 'Data/StJudesOncotreeCodes.RDS'
if(file.exists(tgtFile)){
  url ='http://oncotree.mskcc.org/api/tumorTypes?version=oncotree_latest_stable'
  codes = httr::content(httr::GET(url))
  names(codes) = sapply(codes,function(e) e$code)
  saveRDS(codes,tgtFile)
}else{
  codes = readRDS(tgtFile)
}
#Drop blood cancers and unknowns
mDat = mDat[!sapply(codes[mDat$attr_oncotree_disease_code],function(e) e$tissue) %in% c('Lymphoid','Myeloid','Other'),]
#And drop cancers without much info
mDat = mDat[!mDat$attr_oncotree_disease_code %in% c('ASPS','ES','GIST','IFS','MPNST','SYNS'),]
mDat$disease = sapply(codes[mDat$attr_oncotree_disease_code],function(e) e$name)
#Finalise
rownames(mDat) = mDat$sample_name
judes$exposures = judes$exposures[,mDat$sample_name]
judes$gof = judes$gof[,mDat$sample_name]
judes$raw = judes$raw[,mDat$sample_name]
judes$mDat = mDat
###############
# TCGA cleanup
#Load the TCGA data
dat = list.files(file.path(srcCS,'TCGA','results_TCGA'))
dat = grep('EE_',dat,invert=TRUE,value=TRUE)
#Get the cancer types
cTypes = unique(gsub('_(Fetal|TS|ee)_.*','',dat))
tcga = list()
for(cType in cTypes){
  message(cType)
  #Load early embryo
  tcga[[cType]] = csa$normaliseExposures(file.path(srcCS,'TCGA','results_TCGA',paste0(cType,'_Fetal_new_TS_fitExposures.tsv')))
}
#Now separate out blastoids
blastoid = tcga$HSC_blast
tcga = tcga[-match('HSC_blast',names(tcga))]
#And finalise as they're easy
blastoid$mDat = data.frame(row.names = colnames(blastoid$exposures),
                           sampleID = colnames(blastoid$exposures),
                           disease = 'Blastoid')
#Flatten TCGA, maintaining metadat
mDat = data.frame(sampleID = unlist(lapply(tcga,function(e) colnames(e$exposures))),
                  disease = rep(names(tcga),sapply(tcga,function(e) ncol(e$exposures))))
rownames(mDat) = mDat$sampleID
tcga = list(exposures = do.call(cbind,lapply(tcga,function(e) e$exposures)),
           gof = do.call(cbind,lapply(tcga,function(e) as.matrix(e$gof))),
           raw = do.call(cbind,lapply(tcga,function(e) as.matrix(e$raw))))
#Now cleanup TCGA based on metadata
mani = readRDS('Data/manifestTCGA.RDS')
mani$sampleID = gsub('-','.',mani$UniqueSampleID)
cmnCols = c('dataset','source','tissue','disease','diseaseBroad','biopsy','subtype','outcome','UniqueSampleID','sampleID','tgtFile','toKeep')
mani = mani[match(colnames(tcga$exposures),mani$sampleID),]
#Drop non-primary tumours
mani = mani[mani$biopsy=='Primary Tumor',]
mani$diseaseBroad = mDat$disease[match(mani$sampleID,mDat$sampleID)]
mDat = mani[,cmnCols]
rownames(mDat) = mani$sampleID
tcga$exposures = tcga$exposures[,mDat$sampleID]
tcga$gof = tcga$gof[,mDat$sampleID]
tcga$raw = tcga$raw[,mDat$sampleID]
tcga$mDat = mDat
tcga$bigMani = mani
###########################
# Get reference annotation
dat = list.files(file.path(srcCS,'Single_Cell_Reference_Signatures/references/'),full.names=TRUE)
noms = gsub('\\.tsv$','',basename(dat))
sigs = lapply(dat,read.table,sep='\t',header=TRUE)
sigs = lapply(sigs,colnames)
names(sigs) = noms
dd = data.frame(sig=rownames(tcga$exposures))
rownames(dd) = dd$sig
for(nom in noms)
  dd[,nom] = dd$sig %in% sigs[[nom]]
fetalSigs = dd$sig[dd$cellSigs_Han_wo_Heart_ee_gast & !dd$cellSigs_summary_TS_norm_ee_gast]
adultSigs = dd$sig[!dd$cellSigs_Han_wo_Heart_ee_gast & dd$cellSigs_summary_TS_norm_ee_gast]
eeSigs = dd$sig[dd$earlyemb_summary]
gastSigs = dd$sig[!dd$earlyemb_summary & dd$ee_gastr_summary]
scores = data.frame(sampleID = c(colnames(target$gof),colnames(judes$gof),colnames(tcga$gof),colnames(blastoid$gof)),
                    dataset = rep(c('TARGET','StJudes','TCGA','Blastoid'),c(ncol(target$gof),ncol(judes$gof),ncol(tcga$gof),ncol(blastoid$gof))),
                    disease = c(target$mDat$disease,judes$mDat$disease,tcga$bigMani$cgc_file_disease_type,blastoid$mDat$disease),#The TCGA is the version Maria  used, more detailed is tcga$mDat$disease
                    fetalScore = c(colSums(target$exposures[fetalSigs,]),
                                   colSums(judes$exposures[fetalSigs,]),
                                   colSums(tcga$exposures[fetalSigs,]),
                                   colSums(blastoid$exposures[fetalSigs,])),
                    adultScore = c(colSums(target$exposures[adultSigs,]),
                                   colSums(judes$exposures[adultSigs,]),
                                   colSums(tcga$exposures[adultSigs,]),
                                   colSums(blastoid$exposures[adultSigs,])),
                    eeScore = c(colSums(target$exposures[eeSigs,]),
                                colSums(judes$exposures[eeSigs,]),
                                colSums(tcga$exposures[eeSigs,]),
                                colSums(blastoid$exposures[eeSigs,])),
                    gastScore = c(colSums(target$exposures[gastSigs,]),
                                  colSums(judes$exposures[gastSigs,]),
                                  colSums(tcga$exposures[gastSigs,]),
                                  colSums(blastoid$exposures[gastSigs,])),
                    intercept = c(as.numeric(target$exposures['Intercept',]),
                                  as.numeric(judes$exposures['Intercept',]),
                                  as.numeric(tcga$exposures['Intercept',]),
                                  as.numeric(blastoid$exposures['Intercept',])),
                    pR2 = c(as.numeric(target$gof['pR2',]),
                                  as.numeric(judes$gof['pR2',]),
                                  as.numeric(tcga$gof['pR2',]),
                                  as.numeric(blastoid$gof['pR2',])))
#Make the boxplots in the style of this paper https://www.nature.com/articles/nature25480
dd = scores[scores$dataset!='Blastoid',]
dd = scores
dd$cancerType = factor(ifelse(dd$dataset=='TCGA','Adult',ifelse(dd$dataset=='Blastoid','Control','Childhood')),levels=c('Childhood','Adult','Control'))
#Make disease names unique
dd$disease[dd$disease=='Neuroblastoma' & dd$dataset=='TARGET'] = 'Neurblastoma (TARGET)'
dd$disease[dd$disease=="Wilms Tumor"]  = 'Wilms Tumor (TARGET)'
dd$disease[dd$disease=="Rhabdoid Tumor"]  = 'Rhabdoid Tumor (TARGET)'
#Fix naming for consistency with Maria's
w = dd$dataset=='StJudes'
dd$disease[w] = c(`Adrenocortical Carcinoma` = 'Adrenocortical Tumor',
               `Wilms' Tumor` = 'Wilms Tumor',
               `CNS/Brain` = 'CNS/Brain Tumor',
               `Desmoplastic Small-Round-Cell Tumor` = 'Desmoplastic small round cell Tumor',
               `Ependymomal Tumor` = 'Ependymoma',
               `Choroid Plexus Carcinoma` = 'Choroid Plexus Carcinoma',
               `High-Grade Glioma, NOS` = 'High Grade Glioma',
               `Low-Grade Glioma, NOS` = 'Low Grade Glioma',
               `Medulloblastoma` = 'Medulloblastoma',
               `Melanoma` = 'Melanoma',
               `Neuroblastoma` = 'Neuroblastoma',
               `Osteosarcoma` = 'Osteosarcoma',
               `Retinoblastoma` = 'Retinoblastoma',
               `Rhabdomyosarcoma` = 'Rhabdomyosarcoma')[dd$disease[w]]
#Drop ones that aren't meaningfully "childhood"
dd = dd[!(dd$cancerType=='Childhood' & dd$disease %in% c('Melanoma','Chordoma')),]
####Fix up disease codes
###dd$diseaseShort = dd$disease
###dd$diseaseShort[dd$dataset=='TARGET'] = c(Neuroblastoma='TARGET.NBL',`Rhabdoid Tumor`='TARGET.RT',`Wilms Tumor`='TARGET.WT')[dd$diseaseShort[dd$dataset=='TARGET']]
###dd$diseaseShort[dd$dataset=='StJudes'] = c(`Rhabdomyosarcoma` = 'RHB',`Adrenocortical Tumor` = 'ACT',`Adrenocortical Carcinoma` = 'ACT',`Desmoplastic small round cell Tumor` = 'DSRCT',`Desmoplastic Small-Round-Cell Tumor`='DSRCT',`Osteosarcoma` = 'OS',`Melanoma` = 'MEL',`Chordoma` = 'ST',`Wilms Tumor` = 'SJ.WLM',`Wilms' Tumor`='SJ.WLM',`Neuroblastomas` = 'NBL',`NeuroblastomaSJ`='SJ.NBL',`Medulloblastoma` = 'MB',`High Grade Glioma` = 'HGG',`High-Grade Glioma, NOS`='HGG',`Brain Tumor` = 'BT',`CNS/Brain`='BT',`Ependymoma` = 'EPD',`Ependymomal Tumor`='EPD',`Low Grade Glioma` = 'LGG',`Low-Grade Glioma, NOS`='LGG',`Retinoblastoma` = 'RB',`Choroid Plexus Carcinoma` = 'CPC')[dd$diseaseShort[dd$dataset=='StJudes']]
#Split by disease
ss = split(dd$fetalScore,dd$disease)
#Order by type, then median score 
medScore = sapply(split(dd$fetalScore,dd$disease),median)
cancType = dd$cancerType[match(names(ss),dd$disease)]
ss = ss[order(cancType,medScore)]
######################
# It's plotting time!
pdf(file.path(plotDir,'Fetalness.pdf'),width=7,height=4)
par(mar=c(9,5,1,1))
plot(NA,
     xlim=c(0.5,length(ss)+.5),
     ylim=c(0,1),
     xlab='',
     ylab='Fetalness',
     xaxt='n',
     type='n')
##########
# Shading
cols=list(Childhood='#75335C',
       Adult='#869E48',
       Control='#FC8E62',
       pts='#000000')
cexs = list(pts=0.2)
lwds = list(median=2,
            iqr=0.5)
alphas = list(pts=0.4)
tgtPts=20
boxQuants = c(0.25,0.5,0.75)
left=0.5
for(tgt in levels(cancType)){
  right = left+sum(cancType==tgt)
  xRange=c(left,right)
  quants = quantile(dd$fetalScore[dd$cancerType==tgt],boxQuants)
  rect(xleft=xRange[1],
       xright=xRange[2],
       ybottom=quants[2]-1.5*diff(quants[c(1,3)]),
       ytop=quants[2]+1.5*diff(quants[c(1,3)]),
       border=NA,
       col=paste0(cols[[tgt]],as.hexmode(round(255*0.07))))
  #Q1 - Q3
  rect(xleft=xRange[1],
       xright=xRange[2],
       ybottom=quants[1],
       ytop=quants[3],
       border=NA,
       col=paste0(cols[[tgt]],as.hexmode(round(255*0.3))))
  lines(xRange,rep(quants[2],2),col=cols[[tgt]],lwd=lwds$median)
  left = right
}
###########
# Boxplots
#First points with jitter and target alpha
for(i in seq_along(ss)){
  nom = names(ss)[i]
  dat = ss[[i]]
  #Target this formula being true tgtPts = alpha*nPts
  alpha = as.character(as.hexmode(round(alphas$pts*min(1,tgtPts/length(dat))*255)))
  if(nchar(alpha)==1)
    alpha = paste0('0',alpha)
  points(i+rnorm(length(dat),sd=0.2),dat,col=paste0(cols$pts,alpha),pch=19,cex=cexs$pts)
  lines(i+c(-0.4,0.4),rep(median(dat),2),
        col='black',
        lwd=lwds$median2)
  lines(rep(i,2),quantile(dat,c(0.25,0.75)),
        col='black',
        lwd=lwds$iqr)
}
#Add the axis label ticks
axis(1,
     at=seq_along(ss),
     labels=names(ss),
     las=2,
     gap.axis=-1,
     cex.axis=0.8)
#And the labels themselves
#text(x = seq_along(ss),
#     y = -.08,
#     adj=1.0,
#     labels = names(ss),
#     xpd = NA,
#     ## Rotate the labels by 35 degrees.
#     srt = 35,
#     cex = .2)
dev.off()

#############################
# Bulk transcriptomes table #
#############################

tab = dd[,c('sampleID','dataset','disease','cancerType')]
############################
# Get the extra liver bulks
srcLiverCS = 'Results/cellSignalOutput/Gerdas/liverFull'
liver_meta = read.table(file.path(srcLiverCS,'liver_bulk_sample_metadata.txt'),sep='\t',header=TRUE)
rownames(liver_meta)=liver_meta$sample_name_no_dashes
liver_meta_rel=liver_meta[colnames(liver_cellsig$exposures),]
#The extra ones are the control (GTEx) and hepatoblastomas
t1 = liver_meta_rel[liver_meta_rel$dataset %in% c('GTEx (brain)','Sekiguchi'),]
t1$sampleID = t1$sample_name
t1$disease = t1$sample_type
t1$cancerType = ifelse(t1$dataset=='Sekiguchi','Childhood','Control')
t1 = t1[,colnames(tab)]
tab = rbind(tab,t1)
#Add in fetal liver samples (all 2 of them)
t1 = data.frame(sampleID=c('fetal_Liver_1','fetal_Liver_2'),
                dataset='Gerrard',
                disease='NormalFetalLiver',
                cancerType='Control')
tab = rbind(tab,t1)
#################################
# Get the extra colorectal bulks
srcGutCS = 'Results/cellSignalOutput/Gerdas/gutFull'
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
paths$labs = c(CANCER='Colorectal cancer',
               'NORMAL EPITH' = 'Normal adult gut',
               'Solid Tissue Normal' = 'Normal adult gut',
               'VILLOUS ADENOMA' = 'Villous adenoma',
               'fetal_bulk' = 'Normal fetal gut',
               'gtex_brain' = 'GTEX brain',
               'Primary Tumor' = 'Colorectal cancer')[paths$sample_type]
paths$sampleID = gsub('\\.tsv$','',paths$sample)
paths$disease = paths$labs
paths$dataset[paths$dataset=='fetal bulk']='Kraiczy'
paths$cancerType = ifelse(paths$dataset=='Kraiczy','Control','Adult')
t1 = paths[,colnames(tab)]
t1 = t1[t1$dataset %in% c('Kraiczy','Boardman'),]
t1 = t1[!is.na(t1$disease),]
tab = rbind(tab,t1)
#Fix up references
#NOTE:FILL IN MISSING HERE.
tab$source = tab$dataset
tab$source[tab$source=='Blastoid'] = 'Kagawa et al. 2022'
tab$source[tab$source=='Boardman'] = 'Druliner et al. 2018'
tab$source[tab$source=='Gerrard'] = 'Gerrard et al. 2016'
tab$source[tab$source=='GTEx (brain)'] = 'GTEx Consortium et al. 2017'
tab$source[tab$source=='Kraiczy'] = 'Kraiczy et al. 2019'
tab$source[tab$source=='Sekiguchi'] = 'Sekiguchi et al. 2020'
tab$source[tab$source=='StJudes'] = 'McLeod et al. 2021'
tab$source[tab$source=='TARGET'] = 'https://ocg.cancer.gov/programs/target/data-matrix'
tab$source[tab$source=='TCGA'] = 'https://www.cancer.gov/tcga'
###############
# Write it out
write.table(tab,file.path(plotDir,'bulkTranscriptomes.tsv'),row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

####################################
# Single cell transcriptomes table #
####################################

#Want something like name, description, cell type number, source
tab = list(earlyEmbryo = list(description = 'Cells from human pre-implantation embryo',
                              cellTypes = 2,
                              source = 'Mole et al. 2021'),
           gastrulation = list(description = 'Human cells from three germ layers at gastrulatian',
                               cellTypes = 3,
                               source = 'Tyser et al. 2021'),
           panFetal = list(description = 'Human cells from multiple organs during development',
                           cellTypes = 587,
                           source = 'Han et al. 2020'),
           panAdult = list(description = 'Human cells from multiple post-natal organs',
                           cellTypes = 328,
                           source = 'The Tabula Sapiens Consortium, 2022'),
           fetalLiver = list(description = 'Detailed cell atlas of the human fetal liver',
                             cellTypes = 20,
                             source = 'Popescu et al. 2019'),
           adultLiver = list(description = 'Detailed cell atlas of the human post-natal liver',
                             cellTypes = 13,
                             source = 'MacParland et al. 2018'),
           gut = list(description = 'Detailed cell atlas of the human fetal and post-natal gastrointestional tract.',
                      cellTypes = 33,
                      source = 'Elmentaite et al. 2021'),
           hepatocellularcarcinoma = list(description='Cells from 8 human HCCs',
                                          cellTypes=NA,
                                          source = 'Ho et al. 2021'),
           hepatoblastoma = list(description='Cells from one hepatoblastoma from [REF] year old boy.',
                                 cellTypes=NA,
                                 source='EGAID_PENDING'),
           colorectalcancer = list(description='Cells from 20 human CRCs',
                                   cellTypes=NA,
                                   source = 'Lee et al. 2020')
)
t1 = data.frame(name = names(tab))
for(nom in names(tab[[1]]))
  t1[,nom] = sapply(tab,function(e) e[[nom]])
tab = t1
write.table(tab,file.path(plotDir,'singleCellTranscriptomes.tsv'),row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

#########################
# Stemness correlation  #
#########################

pdf(file.path(plotDir,'CorrelationStemnessAndEE.pdf'),width=4,height=3)
#Load stemness scores
stemness = read.table('Results/eeAndGastScores.tsv',sep='\t',header=TRUE)
stemness$sampleID = gsub('-','.',stemness$sampleName)
d = stemness
plot(d$epiblastAndHypoblastScoreWithEarlyEmbryoRef,d$stemnessScore,
     pch=19,
     cex=0.1,
     xlab='Early embryo signal',
     ylab='Stemness score',
     col='#00000033')
dev.off()

########################
# Normal cell stemness #
########################

nSamps=100000
stemness = read.table('Data/stemnessGeneWeights.tsv',sep='\t',header=FALSE)
stemness = setNames(stemness[,2],stemness[,1])
#Test stemness score for enrichment for cell cycle genes
gsets = lapply(cc.genes.updated.2019,function(e) intersect(e,names(stemness)))
#Get the mean distribution for randomly sampled genes
null = lapply(gsets,function(e) sapply(seq(nSamps),function(ee) mean(sample(stemness,length(e)))))
#p-values
pVals = sapply(names(gsets),function(e) mean(null[[e]]>=mean(stemness[gsets[[e]]])))
#And plots
pdf(file.path(plotDir,'stemnessProlifPval.pdf'),width=4,height=6)
par(mfrow=c(2,1))
for(nom in names(gsets)){
  hist(null[[nom]],
       xlab='Mean coefficient',
       main=nom,
       n=100)
  abline(v=mean(stemness[gsets[[nom]]]))
}
dev.off()

#Liver first
liver = readRDS('Results/liver_and_cortex.rds')
liver = CellCycleScoring(liver,g2m.features = cc.genes.updated.2019$g2m.genes,s.features = cc.genes.updated.2019$s.genes)
noms = unique(liver@meta.data$annot)
tpms = list()
for(nom in noms){
  message(nom)
  w = which(liver@meta.data$annot==nom)
  phases = liver@meta.data$Phase[w]
  for(phase in unique(phases)){
    ww = which(liver@meta.data$annot==nom & liver@meta.data$Phase==phase)
    t1 = rowSums(liver@assays$RNA@counts[,ww,drop=FALSE])
    tpms[[paste0(nom,'_',phase)]] = t1/sum(t1)*1e6
  }
}
tpms = do.call(cbind,tpms)
sscore = calcStemness(tpms,stemness,convertGeneNames=FALSE)
t1 = split(sscore,gsub('.*_','',names(sscore)))
t1 = t1[c('G1','S','G2M')]
#Fix names version
t2 = lapply(t1,function(e) setNames(e, gsub('_(G2M|S|G1)$','',names(e))))
t3 = Reduce(union,lapply(t2,names))
t4 = data.frame(do.call(cbind,lapply(t2,function(e) e[t3])))
t5 = t4[,2:3]-t4[,1]
pdf(file.path(plotDir,'cellCycleStemnessRelative.pdf'),width=4,height=3)
plot(NA,
     xaxt='n',
     type='n',
     las=2,
     xlim=c(0,2),
     ylim=c(0,0.25),
     xlab='Phase of cell cycle',
     ylab='Difference relative to G1'
     )
gap=0.2
left=gap/2
for(nom in colnames(t5)){
  t2 = t5[,nom]
  right = left+1-gap
  xlab = rnorm(length(t2),mean=.5*(left+right),sd=(1-gap)/2*.25)
  xlab[xlab<left]=left
  xlab[xlab>right]=right
  points(xlab,t2,pch=19,cex=0.1)
  #Median
  lines(c(left,right),rep(median(t2,na.rm=TRUE),2))
  #Quartiles
  lines(rep(.5*(left+right),2),quantile(t2,c(.25,.75),na.rm=TRUE))
  left = left+1
}
axis(1,seq(2)-0.5,colnames(t5))
dev.off()
pdf(file.path(plotDir,'cellCycleStemness.pdf'),width=4,height=3)
plot(NA,
     xaxt='n',
     type='n',
     las=2,
     xlim=c(0,3),
     ylim=c(0,1),
     xlab='Phase of cell cycle',
     ylab='Stemness score'
     )
gap=0.2
left=gap/2
for(nom in names(t1)){
  t2 = t1[[nom]]
  right=left+1-gap
  xlab = rnorm(length(t2),mean=.5*(left+right),sd=(1-gap)/2*.25)
  xlab[xlab<left]=left
  xlab[xlab>right]=right
  points(xlab,t2,
         pch=19,
         cex=0.1)
  #Median
  lines(c(left,right),rep(median(t2),2))
  #Quartiles
  lines(rep(.5*(left+right),2),quantile(t2,c(.25,.75)))
  left = left+1
}
axis(1,seq(3)-0.5,names(t1))
dev.off()




                      








