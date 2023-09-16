source('matts-bits-of-code-master/cellSignalAnalysisV2/cellSignalAnalysis.R')
library(reshape2)
library(ggplot2)
library(ggbeeswarm)
library(viridis)
gut_cellsig=normaliseExposures('/home/jovyan/Dediff/cellsig_gut/out_combined_cortex/OutRun_fitExposures.tsv')
paths=read.table('/home/jovyan/Dediff/cellsig_gut/all_paths.txt', header = F)
paths$dataset=c(rep("TCGA", 715), rep("Boardman", 79), rep("fetal bulk", 8), rep("gtex_brain", 370))
paths$sample=basename(paths$V1)

boardman_meta=read.table('boardman_data/boardman_meta.txt', sep = ",", header = T)
boardman_meta_new=data.frame(sample=colnames(gut_cellsig$exposures)[716:794])
which(boardman_meta_new$sample%in%boardman_meta$sample)
boardman_meta_new$sample_type="x"
boardman_meta_new$sample_type=boardman_meta$sample_type[match(boardman_meta_new$sample, boardman_meta$sample)]
boardman_meta_new$sample_type[which(is.na(boardman_meta_new$sample_type))]="unknown"
boardman_meta_new$sample_type_wide=boardman_meta_new$sample_type
boardman_meta_new$sample_type_wide[grep("VILLOUS", boardman_meta_new$sample_type_wide)]= "VILLOUS ADENOMA"

mDat_col=readRDS('mDat_gut.rds')

tcga_meta=data.frame(sample=colnames(gut_cellsig$exposures)[1:715],
                     sample_type=c(as.character(mDat_col$gdc_cases.samples.sample_type)),
                     stage=c(mDat_col$xml_stage_event_pathologic_stage))

paths$sample_type="x"
paths$sample_type[1:715]=tcga_meta$sample_type
paths$sample_type[716:794]=boardman_meta_new$sample_type_wide
paths$sample_type[795:802]="fetal_bulk"
paths$sample_type[803:nrow(paths)]="gtex_brain"

plotExposures(gut_cellsig, column_split=paths$sample_type, show_column_names=F, show_column_dend=F,
              column_title_rot=90, cluster_column_slices=F,
              height = nrow(gut_cellsig$exposures)*unit(3, "mm"), column_gap=unit(3, "mm"))
master_df=as.data.frame(t(gut_cellsig$exposures))
master_df$sample=rownames(master_df)
master_df$total_adult_signal=rowSums(master_df[,grep("Adult", colnames(master_df))])
master_df$total_fetal_signal=rowSums(master_df[,grep("Fetal", colnames(master_df))])
master_df$Adult_stem_ta=(master_df$`Adult_Stem cells`+master_df$Adult_TA)

master_rel_cols=master_df[,c("sample", "Intercept", "total_fetal_signal", "total_adult_signal", "Adult_stem_ta", "Cortex")]

master_rel_cols$sample_type=paths$sample_type
master_rel_cols$dataset=paths$dataset

table(master_rel_cols$sample_type)
master_rel_cols$sample_type[which(master_rel_cols$sample_type=="CANCER")]="Primary Tumor"
master_rel_cols$sample_type[which(master_rel_cols$sample_type=="NORMAL EPITH")]="Solid Tissue Normal"
master_rel_cols$sample_type[which(master_rel_cols$sample_type=="VILLOUS ADENOMA")]="Villous Adenoma"
master_rel_cols=master_rel_cols[master_rel_cols$sample_type%in%c("fetal_bulk","gtex_brain",
                                                 "Primary Tumor", "Solid Tissue Normal", "Villous Adenoma" ),]
master_melted=melt(master_rel_cols, id.vars = c("sample", "sample_type", "dataset"))

master_melted$sample_type=factor(master_melted$sample_type, levels = c("gtex_brain","fetal_bulk", "Solid Tissue Normal",
                                                                       "Villous Adenoma","Primary Tumor"))

master_melted$dataset[master_melted$dataset=="fetal bulk"]="Kraiczy (Fetal bulk)"
master_melted$dataset[master_melted$dataset=="gtex_brain"]="GTEx (Brain)"
master_melted$dataset=factor(master_melted$dataset, levels=c("Kraiczy (Fetal bulk)","Boardman","TCGA", "GTEx (Brain)"))


ggplot(data = master_melted, aes(x=sample_type,y=value, group=dataset)) +
  geom_quasirandom(mapping = aes(x=sample_type,y=value, group=dataset), dodge.width=.8, shape=19, cex=0.5, alpha=0.2) +
  facet_wrap(. ~variable, ncol = 5) +
  stat_summary(mapping = aes(x = sample_type, y = value, shape=dataset, color=dataset),
               geom = "pointrange",
               fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median, 
               position=position_dodge(width=0.8), cex=.8) + scale_color_manual(values = c("orange", "blue","green", "red")) +  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x="Sample type", y = "Proportion of signal explained in sample")

gut_meta=read.table('gut_bulk_sample_metadata.txt', sep = "\t", header = T)
rownames(gut_meta)=gut_meta$sample_name_no_dashes
gut_cellsig=normaliseExposures('/home/jovyan/Dediff/cellsig_everything/out_gut/OutRun_fitExposures.tsv')

mDat_gut=readRDS('mDat_gut.rds')
gut_cellsig_df=as.data.frame(t(gut_cellsig$exposures))
gut_cellsig_df$best_lab=colnames(gut_cellsig_df)[apply(gut_cellsig_df,1,which.max)]
gut_cellsig_tcga=gut_cellsig_df[1:715,]
gut_cellsig_tcga$sample_name=rownames(gut_cellsig_tcga)
gut_cellsig_tcga$tcga_meta_sample_name=rownames(mDat_gut)
gut_cellsig_tcga$sample_type=mDat_gut$cgc_sample_sample_type

just_tums=mDat_gut

rownames(just_tums)=rownames(gut_cellsig_tcga)
gut_cellsig_tum=gut_cellsig_tcga[which(gut_cellsig_tcga$sample_type%in%c( "Primary Tumor","Recurrent Tumor")),]
tcga_gut_tum_meta=just_tums[rownames(gut_cellsig_tum),]
gut_meta_df=as.data.frame(tcga_gut_tum_meta)

rel_gut_meta=data.frame(sample=rownames(gut_cellsig_tum), 
                        stage=tcga_gut_tum_meta$xml_stage_event_pathologic_stage, 
                        location=tcga_gut_tum_meta$tissue,
                        vital_status=tcga_gut_tum_meta$cgc_case_vital_status, 
                        disease_type=tcga_gut_tum_meta$cgc_file_disease_type,
                        percent_tum_nuc=tcga_gut_tum_meta$cgc_slide_percent_tumor_nuclei,
                        percent_tum_cell=tcga_gut_tum_meta$cgc_slide_percent_tumor_cells,
                        percent_stromal=tcga_gut_tum_meta$cgc_slide_percent_stromal_cells,
                        percent_normal=tcga_gut_tum_meta$cgc_slide_percent_normal_cells,
                        percent_necrosis=tcga_gut_tum_meta$cgc_slide_percent_necrosis,
                        pathologic_t=tcga_gut_tum_meta$cgc_case_pathologic_t,
                        pathologic_n=tcga_gut_tum_meta$cgc_case_pathologic_n,
                        gender=tcga_gut_tum_meta$xml_gender,
                        age_at_diagnosis=tcga_gut_tum_meta$xml_age_at_initial_pathologic_diagnosis, 
                        case_tumor_status=tcga_gut_tum_meta$cgc_case_tumor_status,
                        follow_up_tumour_status=tcga_gut_tum_meta$cgc_follow_up_tumor_status,
                        follow_up_vital_status=tcga_gut_tum_meta$cgc_follow_up_vital_status,
                        case_histological_diagnosis=tcga_gut_tum_meta$cgc_case_histological_diagnosis,
                        new_tum_after_treatment=tcga_gut_tum_meta$cgc_follow_up_new_tumor_event_after_initial_treatment,
                        neopl_status=tcga_gut_tum_meta$xml_person_neoplasm_cancer_status, 
                        has_new_tum_events=tcga_gut_tum_meta$xml_has_new_tumor_events_information, 
                        lymphnodes_positive=tcga_gut_tum_meta$xml_number_of_lymphnodes_positive_by_he, 
                        pretr_cea_level=tcga_gut_tum_meta$xml_preoperative_pretreatment_cea_level, 
                        venous_invasion=tcga_gut_tum_meta$xml_venous_invasion, 
                        lymphatic_invasion=tcga_gut_tum_meta$xml_lymphatic_invasion, 
                        perineural_invasion=tcga_gut_tum_meta$xml_perineural_invasion_present, 
                        kras_mutation=tcga_gut_tum_meta$xml_kras_mutation_found, 
                        braf_analysis=tcga_gut_tum_meta$xml_braf_gene_analysis_result, 
                        stem_signal=gut_cellsig_tum$`Adult_Stem cells`, 
                        ta_signal=gut_cellsig_tum$Adult_TA, 
                        intercept=gut_cellsig_tum$Intercept,
                        gast_sig=gut_cellsig_tum$Endoderm+gut_cellsig_tum$Mesoderm+gut_cellsig_tum$`Non-Neural Ectoderm`)

rel_gut_meta$stage[which(rel_gut_meta$stage%in%c("Stage I","Stage IA"))]="Stage I"
rel_gut_meta$stage[which(rel_gut_meta$stage%in%c("Stage IIA","Stage II" ,"Stage IIB","Stage IIC"))]="Stage II"
rel_gut_meta$stage[which(rel_gut_meta$stage%in%c("Stage IIIA","Stage III" ,"Stage IIIB","Stage IIIC"))]="Stage III"
rel_gut_meta$stage[which(rel_gut_meta$stage%in%c("Stage IVA","Stage IV" ,"Stage IVB"))]="Stage IV"

high_na_cols_gut=names(apply(rel_gut_meta, 2, function(x){sum(is.na(x))})[apply(rel_gut_meta, 2, function(x){sum(is.na(x))})>0.2*nrow(rel_gut_meta)])


lm2=lm(stem_signal ~ stage + location + vital_status +disease_type + percent_tum_nuc + percent_tum_cell + percent_stromal +
         percent_normal + percent_necrosis +pathologic_t +pathologic_n + gender + age_at_diagnosis +case_tumor_status +
         case_histological_diagnosis + neopl_status + has_new_tum_events + lymphnodes_positive + venous_invasion +lymphatic_invasion, data = rel_gut_meta)

summary(lm2)
lm3=lm(gast_sig ~ stage + location + vital_status +disease_type + percent_tum_nuc + percent_tum_cell + percent_stromal +
         percent_normal + percent_necrosis +pathologic_t +pathologic_n + gender + age_at_diagnosis +case_tumor_status +
         case_histological_diagnosis + neopl_status + has_new_tum_events + lymphnodes_positive + venous_invasion +lymphatic_invasion, data = rel_gut_meta)
summary(lm3)

ggplot(rel_gut_meta, aes(x=case_histological_diagnosis, y=stem_signal)) +geom_boxplot()

ggplot(rel_gut_meta, aes(x=follow_up_tumour_status, y=gast_sig)) +geom_boxplot()
ggplot(rel_gut_meta, aes(x=new_tum_after_treatment, y=gast_sig)) +geom_boxplot()
ggplot(rel_gut_meta, aes(x=log(pretr_cea_level), y=gast_sig)) +geom_point()
ggplot(rel_gut_meta, aes(x=perineural_invasion, y=gast_sig)) +geom_boxplot()
ggplot(rel_gut_meta, aes(x=kras_mutation, y=gast_sig)) +geom_boxplot()
ggplot(rel_gut_meta, aes(x=braf_analysis, y=gast_sig)) +geom_boxplot()


ggplot(rel_gut_meta, aes(x=follow_up_tumour_status, y=stem_signal)) +geom_boxplot()
ggplot(rel_gut_meta, aes(x=new_tum_after_treatment, y=stem_signal)) +geom_boxplot()
ggplot(rel_gut_meta, aes(x=log(pretr_cea_level), y=stem_signal)) +geom_point()
ggplot(rel_gut_meta, aes(x=perineural_invasion, y=stem_signal)) +geom_boxplot()
ggplot(rel_gut_meta, aes(x=kras_mutation, y=stem_signal)) +geom_boxplot()
ggplot(rel_gut_meta, aes(x=braf_analysis, y=stem_signal)) +geom_boxplot()
ggplot(rel_gut_meta, aes(x=has_new_tum_events, y=stem_signal)) +geom_boxplot()
