
#script to calculate the correlation between early embryob and stemness score and plot the scatterplot
tcga_stjude_stemness=read.table('/lustre/scratch117/casm/team274/gk14/Dediff/eeAndGastScores.tsv', sep = "\t", header = T)
meta_ee_stemness=readRDS('/lustre/scratch117/casm/team274/gk14/Dediff/eeAndGastScoresMetadata.RDS')
summary(lm(tcga_stjude_stemness$stemnessScore~tcga_stjude_stemness$epiblastAndHypoblastScoreWithEarlyEmbryoRef))
summary(lm(tcga_stjude_stemness$stemnessScore~tcga_stjude_stemness$epiblastAndHypoblastScoreWithGastrulationPlusEarlyEmbryoRef))

cor(tcga_stjude_stemness$stemnessScore, tcga_stjude_stemness$epiblastAndHypoblastScoreWithEarlyEmbryoRef, method = "pearson")
cor(tcga_stjude_stemness$stemnessScore, tcga_stjude_stemness$epiblastAndHypoblastScoreWithGastrulationPlusEarlyEmbryoRef)


plot(density(tcga_stjude_stemness$epiblastAndHypoblastScoreWithEarlyEmbryoRef))
qqnorm(tcga_stjude_stemness$epiblastAndHypoblastScoreWithEarlyEmbryoRef)
qqline(tcga_stjude_stemness$epiblastAndHypoblastScoreWithEarlyEmbryoRef)

wilcox.test(tcga_stjude_stemness$epiblastAndHypoblastScoreWithEarlyEmbryoRef,tcga_stjude_stemness$stemnessScore,paired=T) 
wilcox.test(tcga_stjude_stemness$epiblastAndHypoblastScoreWithGastrulationPlusEarlyEmbryoRef,tcga_stjude_stemness$stemnessScore,paired=T) 

cor.test(tcga_stjude_stemness$epiblastAndHypoblastScoreWithEarlyEmbryoRef,tcga_stjude_stemness$stemnessScore)
plot(tcga_stjude_stemness$epiblastAndHypoblastScoreWithEarlyEmbryoRef,tcga_stjude_stemness$stemnessScore) 

plot(tcga_stjude_stemness$epiblastAndHypoblastScoreWithGastrulationPlusEarlyEmbryoRef,tcga_stjude_stemness$stemnessScore) 

ggplot(tcga_stjude_stemness, aes(y=stemnessScore, x=epiblastAndHypoblastScoreWithEarlyEmbryoRef)) +geom_point(shape=19, alpha=0.05, size=1) +
  stat_cor(label.y = 0.8, aes(label =  ..r.label..))+ #this means at .6th unit in the y axis, the r squared and p value will be shown
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Early embryo (epiblast + hypoblast) exposures", y="Stemness score")

#calculate the change in early embyro signal when gastrulation reference is added for TCGA/stJudes
tcga_stjude_stemness$change_in_signal=tcga_stjude_stemness$epiblastAndHypoblastScoreWithEarlyEmbryoRef-tcga_stjude_stemness$epiblastAndHypoblastScoreWithGastrulationPlusEarlyEmbryoRef
tcga_stjude_stemness$p_or_a="adult"
tcga_stjude_stemness$p_or_a[9963:11160]="pediatric"


tcga_stjude_stemness$disease_type=c(meta_ee_stemness$mDatTCGA$disease, meta_ee_stemness$mDatStJudes$attr_diagnosis)

which(is.na(tcga_stjude_stemness$disease_type))

na_rm=tcga_stjude_stemness[-which(is.na(tcga_stjude_stemness$disease_type)),]
empty_rm=na_rm[-which(na_rm$disease_type%in%c("", "Other, specify")),]
empty_rm$disease_type[grep("Ependymoma", empty_rm$disease_type)]= "Ependymoma, NOS"
empty_rm$disease_type[grep("Retinoblastoma", empty_rm$disease_type)]= "Retinoblastoma, NOS"
empty_rm$disease_type[grep("Osteosarcoma", empty_rm$disease_type)]= "Osteosarcoma, NOS"
empty_rm$disease_type[grep("Melanoma", empty_rm$disease_type)]= "Melanoma, NOS"
empty_rm$disease_type[which(empty_rm$disease_type%in%c("B-cell Acute Lymphoblastic Leukemia,ETV6-RUNX1",
                                                       "B-cell Acute Lymphoblastic Leukemia, ETV6-RUNX1"))]="B-cell Acute Lymphoblastic Leukemia, ETV6-RUNX1"

empty_rm$disease_type[which(empty_rm$disease_type%in%c("B-cell Acute Lymphoblastic Leukemia,ERG",
                                                       "B-cell Acute Lymphoblastic Leukemia, ERG" ))]="B-cell Acute Lymphoblastic Leukemia, ERG"

empty_rm$disease_type[which(empty_rm$disease_type%in%c("B-cell Acute Lymphoblastic Leukemia,E2A-PBX1",
                                                       "B-cell Acute Lymphoblastic Leukemia, E2A-PBX1"))]="B-cell Acute Lymphoblastic Leukemia, E2A-PBX1" 

empty_rm$disease_type[which(empty_rm$disease_type%in%c("B-cell Acute Lymphoblastic Leukemia, NOS",
                                                       "B-cell Acute Lymphoblastic Leukemia,NOS"))]="B-cell Acute Lymphoblastic Leukemia, NOS" 


pvals=sapply(unique(empty_rm$disease_type), function(x){
  df=empty_rm[which(empty_rm$disease_type==x),]
  if (nrow(df)==1) {
    return(NA)
  } else {
    w_out=wilcox.test(df$epiblastAndHypoblastScoreWithEarlyEmbryoRef, df$epiblastAndHypoblastScoreWithGastrulationPlusEarlyEmbryoRef)
    return(w_out$p.value)
  }
})
qvals=p.adjust(pvals)


empty_rm$qval=qvals

      
ggplot(empty_rm, aes(x=reorder(disease_type, -change_in_signal), y=change_in_signal, fill=p_or_a)) +
  geom_boxplot() +coord_flip()




