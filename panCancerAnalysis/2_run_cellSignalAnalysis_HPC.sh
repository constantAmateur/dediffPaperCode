#!/bin/bash

###################################################
#### run cellSignalAnalysis  deconvolution ########
###################################################


metadata=/nfs/research/icortes/mkalyva/Deconvolution/metadata_to_run.tsv

# example metadata file 
# BULK_SAMPLES  cellSignature_reference output_name
# Bladder_TCGA.txt cellSigs_summary_TS_sapiens_norm.tsv  Bladder_TS
# blastoid_launch.tsv  cellSigs_summary_TS_sapiens_norm.tsv  HSC_blast_TS
# Bone_Marrow_TCGA.txt cellSigs_summary_TS_sapiens_norm.tsv  Bone_Marrow_TS
# Brain_TCGA.txt   cellSigs_summary_TS_sapiens_norm.tsv  Brain_TCGA_TS
# Breast_TCGA.txt  cellSigs_summary_TS_sapiens_norm.tsv  Breast_TCGA_TS ...


while read j
do

bulk=$(echo ${j} | awk '{print $1}')
reference=$(echo ${j} | awk '{print $2}')
out=$(echo ${j} | awk '{print $3}')
log=${out}.log

bsub -M 40G -J T_${reference}${bulk} -o ~/o_${log} -e ~/e_${log}  "ipython -i -- /nfs/research/icortes/mkalyva/SOFTWARE/cellSignalAnalysis/cellSignalAnalysis.py \
  -b ${bulk} \
  -s ${reference} \
  -w ~/cellSignalAnalysis/geneWeights.tsv \
   ~/results_TCGA/${out}"


done<${metadata}