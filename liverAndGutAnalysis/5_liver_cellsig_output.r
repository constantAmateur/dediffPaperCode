source('matts-bits-of-code-master/cellSignalAnalysisV2/cellSignalAnalysis.R')
library(reshape2)
library(ggplot2)
library(ggbeeswarm)
library(viridis)
#cell signal analysis output
liver_meta=read.table('liver_bulk_sample_metadata.txt', sep = "\t", header = T)
rownames(liver_meta)=liver_meta$sample_name_no_dashes
liver_cellsig=normaliseExposures('/home/jovyan/Dediff/cellsig_everything/out_liver/OutRun_fitExposures.tsv')

liver_meta_rel=liver_meta[colnames(liver_cellsig$exposures),]
plotExposures(liver_cellsig, column_split=liver_meta_rel$dataset_and_sample_type,
              show_column_names=F, show_column_dend=F,column_title_rot=90, cluster_column_slices=F,
              height = nrow(liver_cellsig$exposures)*unit(2, "mm"), column_gap=unit(3, "mm"))

fetal_liver_cellsig=normaliseExposures('/home/jovyan/Dediff/cellsig_everything/out_fetal_liver/OutRun_fitExposures.tsv')
plotExposures(fetal_liver_cellsig)

liver_cellsig_df=as.data.frame(liver_cellsig$exposures)
fetal_liver_cellsig_df=as.data.frame(fetal_liver_cellsig$exposures)
lc_df=as.data.frame(t(cbind(liver_cellsig_df, fetal_liver_cellsig_df[rownames(liver_cellsig_df),])))

lc_df$total_adult_signal=rowSums(lc_df[,grep("Adult", colnames(lc_df))])
lc_df$total_fetal_signal=rowSums(lc_df[,grep("Fetal", colnames(lc_df))])
lc_df_rel=lc_df[,c("total_adult_signal", "total_fetal_signal", "Intercept", "Cortex", "Adult_Hepatocytes", "Fetal_Hepatocyte")]
lc_df_rel$dataset=liver_meta$dataset
lc_df_rel$sample_type=liver_meta$sample_type

lc_df_rel=lc_df_rel[lc_df_rel$sample_type%in%c("Hepatocellular Carcinoma", "Solid Tissue Normal",
                                               "Hepatoblastoma","GTEx (brain)","Hanley (fetal liver)"),]
lc_df_rel$sample=rownames(lc_df_rel)
lc_df_rel_melted =melt(lc_df_rel, id.vars = c("sample", "sample_type", "dataset"))
lc_df_rel_melted$sample_type[which(lc_df_rel_melted$sample_type=="Hanley (fetal liver)")] = "Normal liver (fetal)"
lc_df_rel_melted$sample_type[which(lc_df_rel_melted$sample_type=="Solid Tissue Normal")] = "Normal liver (adult)"
lc_df_rel_melted$sample_type=factor(lc_df_rel_melted$sample_type, levels= c("GTEx (brain)", "Normal liver (fetal)",  "Normal liver (adult)",
                                                                            "Hepatoblastoma", "Hepatocellular Carcinoma"))

lc_df_rel_melted$variable=as.character(lc_df_rel_melted$variable)
lc_df_rel_melted$variable[which(lc_df_rel_melted$variable=="Intercept")] = "Unexplained signal"
lc_df_rel_melted$variable[which(lc_df_rel_melted$variable=="Cortex")] = "Negative control\n(adrenal cortex)"
lc_df_rel_melted$variable[which(lc_df_rel_melted$variable=="total_fetal_signal")] = "Total\nfetal signal"
lc_df_rel_melted$variable[which(lc_df_rel_melted$variable=="total_adult_signal")] = "Total\nadult signal"
lc_df_rel_melted$variable[which(lc_df_rel_melted$variable=="Fetal_Hepatocyte")] = "Fetal\nhepatocyte"
lc_df_rel_melted$variable[which(lc_df_rel_melted$variable=="Adult_Hepatocytes")] = "Adult\nhepatocyte"
lc_df_rel_melted$variable[which(lc_df_rel_melted$variable=="Unexplained signal")] = "Unexplained\nsignal"
lc_df_rel_melted$variable=factor(lc_df_rel_melted$variable, levels = c("Total\nfetal signal", "Total\nadult signal", 
                                                                       "Fetal\nhepatocyte", "Adult\nhepatocyte", "Unexplained\nsignal",
                                                                       "Negative control\n(adrenal cortex)"))
ggplot(data = lc_df_rel_melted, aes(x=sample_type,y=value)) +
  geom_quasirandom(mapping = aes(x=sample_type,y=value), dodge.width=.8, shape=19, cex=0.5, alpha=0.2) +
  facet_wrap(. ~variable, ncol = 6) +
  stat_summary(mapping = aes(x = sample_type, y = value),
               geom = "pointrange",
               fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median, cex=.8, shape=15, color="#001c55") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.position = "none") + 
  labs(x="Sample type", y = "Contribution to the sample") 



ggplot(data = lc_df_rel_melted, aes(x=sample_type,y=value)) +
  geom_quasirandom(mapping = aes(x=sample_type,y=value), dodge.width=.8, shape=19, cex=0.5, alpha=0.2) +
  facet_wrap(. ~variable, ncol = 4) +
  stat_summary(mapping = aes(x = sample_type, y = value, shape=sample_type),
               geom = "pointrange",
               fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median, 
               position=position_dodge(width=0.8), cex=.8) +
  scale_color_manual(values = c("orange", "blue","green", "red", "black", "purple")) + 
  scale_shape_manual(name = "Group", values = c(15, 15, 15, 15, 15, 15)) +
  labs(x="Sample type", y = "Proportion of signal explained in sample") + theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.position = "none")
  


mDat_liver=readRDS('mDat_liver.rds')
liver_cellsig_df=as.data.frame(t(liver_cellsig$exposures))
liver_cellsig_df$best_lab=colnames(liver_cellsig_df)[apply(liver_cellsig_df,1,which.max)]
liver_cellsig_tcga=liver_cellsig_df[1:424,]
liver_cellsig_tcga$sample_name=rownames(liver_cellsig_tcga)
liver_cellsig_tcga$tcga_meta_sample_name=rownames(mDat_liver)
liver_cellsig_tcga$sample_type=mDat_liver$cgc_sample_sample_type
#get the bulk liver counts and convert to TPMs, then add to the liver_cellsig_TCGA file
rc2 = readRDS('/lustre/scratch119/realdata/mdt1/team274/my4/bulkRNAseq/rse_gene_TCGA.RDS')
liver_se=rc2[,rownames(mDat_liver)]
gene_info=data.frame( length=liver_se@rowRanges$bp_length, kb_length=liver_se@rowRanges$bp_length/1000)
count_mtx=assays(liver_se)$counts
count_mtx_rpk=apply(count_mtx, 2, function(x){x/gene_info$kb_length})
scale_factor=colSums(count_mtx_rpk)/1000000
count_mtx_tpm=t(apply(count_mtx_rpk, 1, function(x){x/scale_factor}))
grep("ENSG00000167701",rownames(count_mtx_tpm))
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

ggplot(rel_meta, aes(y=fetal_hep, x=log(AFP_tpm))) +geom_point(shape=19, alpha=0.2) +
  geom_smooth(method='lm', formula= y~x, color="#001c55") +
  stat_cor(label.y = 0.6, aes(label =  ..rr.label..))+ #this means at 35th unit in the y axis, the r squared and p value will be shown
  stat_regline_equation(label.y = 0.5) + #this means at 30th unit regresion line equation will be shown
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="AFP expression (log(TPM))", y="Fetal hepatocyte signal")
summary(lm(log(AFP_tpm+1)~fetal_hep, data=rel_meta))

rel_meta$stage[rel_meta$stage%in%c("Stage III",  "Stage IIIA", "Stage IIIB", "Stage IIIC")]="Stage III"
rel_meta$stage[rel_meta$stage%in%c("Stage IV",  "Stage IVA", "Stage IVB")]="Stage IV"
rel_meta$histologic_grade=as.character(rel_meta$histologic_grade)
rel_meta$histologic_grade[which(rel_meta$histologic_grade=="G1")]= 1
rel_meta$histologic_grade[which(rel_meta$histologic_grade=="G2")]=2
rel_meta$histologic_grade[which(rel_meta$histologic_grade=="G3")]=3
rel_meta$histologic_grade[which(rel_meta$histologic_grade=="G4")]=4
rel_meta$histologic_grade=as.numeric(rel_meta$histologic_grade)
#code tumor grades as cont. var

rel_meta$rf_alcohol_use="No"
rel_meta$rf_alcohol_use[grepl("Alcohol", rel_meta$risk_factors) & !grepl("Non-Alcoholic", rel_meta$risk_factors)]= "Yes"

rel_meta$rf_hepc="No"
rel_meta$rf_hepc[grep("Hepatitis C", rel_meta$risk_factors)]="Yes"


rel_meta$rf_yesno="Yes"
rel_meta$rf_yesno[grep("No History of Primary Risk Factors", rel_meta$risk_factors)]="No"

high_na_cols=names(apply(rel_meta, 2, function(x){sum(is.na(x))})[apply(rel_meta, 2, function(x){sum(is.na(x))})>0.2*nrow(rel_meta)])


table(rel_meta$family_cancer_history)
colnames(rel_meta)[which(!colnames(rel_meta)%in%high_na_cols)]
lm1=lm(fetal_hep~log(AFP_tpm+1)+vital_status+percent_tum_nuc+
         percent_tum_cell+percent_stromal+percent_normal+
         percent_necrosis+pathologic_t+pathologic_n +
         gender+age_at_diagnosis+case_tumor_status+
         neopl_status+has_new_tum_events+
         ethnicity+bmi+histologic_grade+rf_yesno+
         family_cancer_history, data = rel_meta)

summary(lm1)

ggplot(rel_meta, aes(x=ethnicity, y=fetal_hep)) +geom_boxplot()

ggplot(rel_meta, aes(x=follow_up_vital_status, y=log(AFP_tpm))) +geom_boxplot()

lm2=lm(gast_sig~log(AFP_tpm+1)+vital_status+percent_tum_nuc+
         percent_tum_cell+percent_stromal+percent_normal+
         percent_necrosis+pathologic_t+pathologic_n +
         gender+age_at_diagnosis+case_tumor_status+
         neopl_status+has_new_tum_events+
         ethnicity+bmi+histologic_grade+rf_yesno+
         family_cancer_history, data = rel_meta)

summary(lm2)



ggplot(rel_meta, aes(x=days_to_death, y=fetal_hep)) +geom_point()
ggplot(rel_meta, aes(x=histologic_grade, y=gast_sig)) +geom_boxplot()
ggplot(rel_meta, aes(x=follow_up_vital_status, y=fetal_hep)) +geom_boxplot()
ggplot(rel_meta, aes(x=new_tum_after_treatment, y=fetal_hep)) +geom_boxplot()
ggplot(rel_meta, aes(x=adjacent_inflamation, y=fetal_hep)) +geom_boxplot()
ggplot(rel_meta, aes(x=viral_serologies, y=fetal_hep)) +geom_boxplot()
ggplot(rel_meta, aes(x=fibrosis_score, y=fetal_hep)) +geom_boxplot()

ggplot(rel_meta, aes(x=histologic_grade, y=adult_hhyp)) +geom_boxplot()
ggplot(rel_meta, aes(x=gast_sig, y=fetal_hep)) +geom_point()
summary(lm1)

lm2=lm(fetal_hep~risk_factors, data = rel_meta)
summary(lm(fetal_hep~fibrosis_score, data = rel_meta))
summary(lm2)

rel_meta$rf_alcohol_use="No"
rel_meta$rf_alcohol_use[grepl("Alcohol", rel_meta$risk_factors) & !grepl("Non-Alcoholic", rel_meta$risk_factors)]= "Yes"

rel_meta$rf_hepc="No"
rel_meta$rf_hepc[grep("Hepatitis C", rel_meta$risk_factors)]="Yes"

lm3=lm(fetal_hep~rf_hepc*rf_alcohol_use, data = rel_meta)

summary(lm3)

ggplot(rel_meta, aes(x=rf_hepc, y=fetal_hep)) +geom_boxplot()

rel_meta$rf_yesno="Yes"
rel_meta$rf_yesno[grep("No History of Primary Risk Factors", rel_meta$risk_factors)]="No"

lm3=lm(fetal_hep~rf_yesno, data = rel_meta)
summary(lm3)

ggplot(rel_meta, aes(x=rf_yesno, y=fetal_hep)) +geom_boxplot()


require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)


test <- multinom(best_lab~stage+vital_status+percent_tum_nuc+
                   percent_tum_cell+percent_stromal+percent_normal+
                   percent_necrosis+pathologic_t+pathologic_n +
                   gender+age_at_diagnosis+case_tumor_status+
                   neopl_status+has_new_tum_events+
                   ethnicity+bmi+histologic_grade+
                   family_cancer_history, data = rel_meta)
summary(test)

z <- summary(test)$coefficients/summary(test)$standard.errors
z

p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

pp <- fitted(test)

lm_dtd=lm(fetal_hep~days_to_death, data = rel_meta)
summary(lm_dtd)


rm_low_afp=rel_meta[which(rel_meta$fetoprotein_outcome_value>20),]

ggplot(rm_low_afp, aes(x=days_to_death, y=fetal_hep)) +geom_point() 
ggplot(rm_low_afp, aes(x=log(fetoprotein_outcome_value), y=fetal_hep)) +geom_point() + geom_smooth(method='lm', formula= y~x)
ggplot(rm_low_afp, aes(x=follow_up_tumour_status, y=fetal_hep)) +geom_boxplot()
ggplot(rm_low_afp, aes(x=follow_up_vital_status, y=fetal_hep)) +geom_boxplot()
ggplot(rm_low_afp, aes(x=new_tum_after_treatment, y=fetal_hep)) +geom_boxplot()
ggplot(rm_low_afp, aes(x=adjacent_inflamation, y=fetal_hep)) +geom_boxplot()
ggplot(rm_low_afp, aes(x=viral_serologies, y=fetal_hep)) +geom_boxplot()
ggplot(rm_low_afp, aes(x=fibrosis_score, y=fetal_hep)) +geom_boxplot() + theme(axis.text.x = element_text(angle = 90))

ggplot(rm_low_afp, aes(x=histologic_grade, y=fetal_hep)) +geom_boxplot()

lm3=lm(fetal_hep~stage+vital_status+percent_tum_nuc+
         percent_tum_cell+percent_stromal+percent_normal+
         percent_necrosis+pathologic_t+pathologic_n +
         gender+age_at_diagnosis+case_tumor_status+
         neopl_status+has_new_tum_events+
         ethnicity+bmi+histologic_grade+
         family_cancer_history +rf_alcohol_use+rf_hepc, data = rel_meta)

summary(lm3)
