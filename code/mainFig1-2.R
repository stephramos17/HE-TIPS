setwd("/mnt/zarrinpar/scratch/sfloresr/HE/HE_metab/liverpath/")

library("tidyverse")
library("data.table")
library("qiime2R")
library("ggpubr")
library("VennDiagram")
library("ape")
library("ggsci")
library("phyloseq")

library(vegan)
# library("RColorBrewer")
# library("DivNet")
#library("vegan")

##################################################
md<-fread("/mnt/zarrinpar/scratch/sfloresr/HE/HE_metab/data_files/metadata_table/metadata_table-00000.tsv")

# md$ATTRIBUTE_blood_procedure[md$ATTRIBUTE_Sample=="Pre-HET"]<-"pre"
# md$ATTRIBUTE_blood_procedure[grepl("Post-",md$ATTRIBUTE_Sample)]<-"post"
# md$ATTRIBUTE_blood_procedure[grepl("POD-1",md$ATTRIBUTE_Sample)]<-"post"
md$ATTRIBUTE_blood_procedure<-factor(md$ATTRIBUTE_blood_procedure,c("pre","post","bd","return"))
md<-md%>%rename(sample_name=filename)%>%rename(sampleid=ATTRIBUTE_sampleID)%>%select(sampleid,everything())%>%select(-sample_name)%>%
  mutate(ATTRIBUTE_blood_procedure_num=ifelse(ATTRIBUTE_blood_procedure=="pre",0,
                                              ifelse(ATTRIBUTE_blood_procedure=="post",1,2)))%>%
  mutate(ATTRIBUTE_pt_id_new=paste(ATTRIBUTE_pt_id,ATTRIBUTE_blood_origin,sep="_"))
write.table(md,"/data/Stephany/liverpath/data_files/LCMS/metadata_table/metadata_table-clean.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("data_files/LCMS/metadata_table/metadata_table-clean.txt")%>%group_by(ATTRIBUTE_pt_id)%>%
  select(ATTRIBUTE_Sex,ATTRIBUTE_Race.Ethnicity,ATTRIBUTE_Age,ATTRIBUTE_diagnosis_alcohol,ATTRIBUTE_diagnosis_nash)%>%unique()%>%
  rename(sampleid=ATTRIBUTE_pt_id)
write.table(md,"/data/Stephany/liverpath/data_files/LCMS/metadata_table/metadata_table-subject.txt",sep = "\t",row.names = FALSE,quote=FALSE)


md<-fread("/data/Stephany/liverpath/data_files/LCMS/metadata_table/metadata_table-clean.txt")%>%
  select(-ATTRIBUTE_blood_procedure_num,-ATTRIBUTE_blood_origin)
logrt<-fread("/data/Stephany/liverpath/diversity_metrics_LCMS/ctf-results/sample_plot_data_procedure.tsv")%>%
  rename(sampleid=`Sample ID`)%>%left_join(.,md,by="sampleid")%>%unique()
write.table(logrt,"/data/Stephany/liverpath/diversity_metrics_LCMS/ctf-results/merged_sample_plot_data_procedure.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("/data/Stephany/liverpath/data_files/LCMS/metadata_table/metadata_table-post.txt")%>%
  select(-ATTRIBUTE_Worst_PostTIPS_HE,-ATTRIBUTE_blood_origin)

logrt<-fread("/data/Stephany/liverpath/diversity_metrics_LCMS/ctf-results-HE/sample_plot_data.tsv")%>%
  rename(sampleid=`Sample ID`)%>%left_join(.,md,by="sampleid")%>%unique()
write.table(logrt,"/data/Stephany/liverpath/diversity_metrics_LCMS/ctf-results-HE/merged_sample_plot_data.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#peripheral-no bd
md<-fread("data_files/LCMS/metadata_table/metadata_table-00000.tsv")

# md$ATTRIBUTE_blood_procedure[md$ATTRIBUTE_Sample=="Pre-HET"]<-"pre"
# md$ATTRIBUTE_blood_procedure[grepl("Post-",md$ATTRIBUTE_Sample)]<-"post"
# md$ATTRIBUTE_blood_procedure[grepl("POD-1",md$ATTRIBUTE_Sample)]<-"post"
md$ATTRIBUTE_blood_procedure<-factor(md$ATTRIBUTE_blood_procedure,c("pre","post","bd","return"))
md<-md%>%dplyr::rename(sample_name="filename")%>%
  filter(ATTRIBUTE_blood_origin=="peripheral",(ATTRIBUTE_blood_procedure=="pre"|ATTRIBUTE_blood_procedure=="post"))%>%dplyr::rename(sampleid=ATTRIBUTE_sampleID)%>%select(sampleid,everything())%>%select(-sample_name)
write.table(md,"data_files/LCMS/metadata_table/metadata_table-peripheral-prepost.txt",sep = "\t",row.names = FALSE,quote=FALSE)

perif_ord <- read_qza("diversity_metrics_LCMS/ordination-peripheral-prepost.qza")

rpca<-perif_ord$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  dplyr::rename(sampleid="SampleID")%>%
  left_join(md,by="sampleid")

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=ATTRIBUTE_blood_procedure)) +
  geom_point(alpha=1.0) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +stat_ellipse(type = "t", linetype = 2,aes(group = ATTRIBUTE_blood_procedure))+
  scale_color_manual(values=c("#52A7DC","#DF643B","#2D653D"))+
  labs(x =paste("PC1 (",round(perif_ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(perif_ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       color="TIPS visit")+ggtitle("Peripheral")+ theme(plot.title = element_text(face = "bold"))

#ggsave("diversity_metrics_LCMS/SFR21_1025_PCoAvisitPerifOnlyPrePostAnnot_wellipse.tiff", plot=p,height=3, width=3)
ggsave("diversity_metrics_LCMS/SFR22_0512_PCoAvisitPerifOnlyPrePostAnnot_wellipse.pdf", plot=p,height=3, width=3)
#ggsave("diversity_metrics_LCMS/SFR21_1006_PCoAvisitPerifOnlyPrePostAnnot_braycurtis_wellipse.tiff", plot=p,height=3, width=3)

mtb_perif<-fread("diff_HE_analysis/SFR21_1025_HEmtbdata_allwmd.txt")%>%
  mutate(feature_name=ifelse(is.na(Compound_Name),paste("mtb ",rowID," m/z ",signif(`row m/z`,digits=4),sep=""),paste("mtb ",rowID," m/z ",signif(`row m/z`,digits=4)," ",Compound_Name,sep="")))%>% 
           dplyr::rename(FeatureID=rowID)%>%
           select(FeatureID, superclass, class, npclassifier_superclass, npclassifier_class,npclassifier_pathway,Compound_Name,feature_name)%>%unique() #%>%column_to_rownames("FeatureID")
subclass<-fread("data_files/LCMS/DB_result/20211105_all_annot_subclass.csv")%>%select(2,1,3,4)
mtb_perif<-left_join(mtb_perif,subclass,by="FeatureID")%>%
   mutate(FeatureID=paste("s",FeatureID,sep=""))
write.table(mtb_perif,"data_files/LCMS/DB_result/spectral_matches_clean_mmvec.txt",sep = "\t",row.names = FALSE,quote=FALSE)

baseplot<-
  ggplot() +
  theme_bw() +
  labs(x =paste("PC1 (",round(perif_ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(perif_ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       fill="TIPS visit")+ggtitle("Peripheral")+
  geom_segment(data=perif_ord$data$Species %>% 
                 mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
                 top_n(5, a) %>% #keep 8 furthest away points
                 mutate(PC1=0.19*PC1, PC2=0.39*PC2)%>% # scale arrows linearly... is this ok? 
                 left_join(mtb_perif,by="FeatureID"),
               aes(x=0, xend=PC1, y=0, yend=PC2,color=feature_name),
               arrow = arrow(length = unit(0.3,"cm"))
  )

# now overlay samples
baseplot<-baseplot +
  geom_point(
    data=rpca,
    aes(x=PC1, y=PC2, fill=ATTRIBUTE_blood_procedure), 
    shape=21
  ) + scale_fill_manual(values=c("#52A7DC","#DF643B","#2D653D")) +theme(plot.title = element_text(face = "bold"))
ggsave("diversity_metrics_LCMS/SFR22_0228_PCoAvisitPerifOnlyPrePostAnnot_warrows.tiff", plot=baseplot,height=4, width=19)

#10% highest/lowest ranked features -natural log
df_div<-fread("diversity_metrics_LCMS/ordination-peripheral-prepost/sample_plot_data_peripheral_prepost.tsv")%>%select(1:3)%>%
  mutate(ATTRIBUTE_blood_procedure=factor(ATTRIBUTE_blood_procedure, levels = c("pre", "post","bd")))

p<-ggplot(df_div, aes(x=ATTRIBUTE_blood_procedure, y=Current_Natural_Log_Ratio,fill=ATTRIBUTE_blood_procedure)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#52A7DC","#DF643B","#2D653D"))+scale_fill_manual(values=c("#52A7DC","#DF643B","#2D653D"))+
  labs(x="TIPS visit",y="Natural Log Ratio")+
  theme(legend.position = "none")

ggsave("diversity_metrics_LCMS/ordination-peripheral-prepost//SFR22_0301_deicode_natlog.pdf",p, height=3
       , width=2.3)

pairwise.wilcox.test(df_div$Current_Natural_Log_Ratio, df_div$ATTRIBUTE_blood_procedure,
                     p.adjust.method="fdr")
# pre  
# post 0.023

rank_feat<-fread("diversity_metrics_LCMS/ordination-peripheral-prepost/selected_features_peripheral_prepost.tsv")
#plot the lowest ranked BAs & AA
rpca<-perif_ord$data$Species %>%
  select(FeatureID, PC1, PC2) %>%
  left_join(mtb_perif,by="FeatureID")%>%
  filter(FeatureID %in% rank_feat$FeatureID)%>%
  mutate(feature_name=fct_reorder(feature_name,PC1))%>%
  filter(annotation=="Annotated")

p<-ggplot(data=rpca, aes(x=feature_name, y= PC1)) +
  geom_bar(stat="identity",fill="#A9A9A9")+ theme_minimal()+coord_flip()+
  labs(x="Feature",y="PC1",title="highest and lowest ranking features")
ggsave("diversity_metrics_LCMS/ordination-peripheral-prepost/SFR22_0304_topbottom10p_PC1.pdf",p, height=8, width=23.5)

HEdf_Perif<-fread("data_files/LCMS/quantification_table/quantification_table-00000.csv")%>%
  select(`row ID`, all_of(md$sampleid))
write.table(HEdf_Perif,"/data/Stephany/liverpath/data_files/LCMS/quantification_table/quantification_table-peripheral-prepost.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#HEdf_dc<-fread("diversity_metrics_LCMS/beta_deicode_lcms_peripheral-prepost/distance-matrix.tsv")%>%column_to_rownames("V1")\
HEdf_dc<-fread("diversity_metrics_LCMS/beta_braycurtis_lcms_peripheral_prepost/distance-matrix.tsv")%>%column_to_rownames("V1")
HRdc_pcoa <- pcoa(HEdf_dc)
HRdc_pcoa_df<-as.data.frame(HRdc_pcoa$vectors)%>%rownames_to_column("sampleid")

p<-HRdc_pcoa_df %>%
  select(sampleid, Axis.1,Axis.2) %>%
  left_join(.,md, by="sampleid") %>%
  ggplot(aes(x=Axis.1, y=Axis.2, color=ATTRIBUTE_blood_procedure)) +
  geom_point(alpha=1.0) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +stat_ellipse(type = "t", linetype = 2,aes(group = ATTRIBUTE_blood_procedure))+
  scale_color_manual(values=c("#52A7DC","#DF643B","#2D653D"))+
  labs(x =paste("PC1 (",round(HRdc_pcoa$values$Relative_eig[1]*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(HRdc_pcoa$values$Relative_eig[2]*100,digits=2),"%)",sep=""),
       color="TIPS visit")+ggtitle("Peripheral") + theme(plot.title = element_text(face = "bold"))
#ggsave("diversity_metrics_LCMS/SFR21_0824_PCoAvisitPerifOnlyPrePostAnnot_wellipse.tiff", plot=p,height=3, width=3)
ggsave("diversity_metrics_LCMS/SFR21_1006_PCoAvisitPerifOnlyPrePostAnnot_braycurtis_wellipse.tiff", plot=p,height=3, width=3)

adonis2(formula = HEdf_dc ~ ATTRIBUTE_blood_procedure, data = md, permutations = 999, method = "braycurtis", by = NULL)

set.seed(1234)
# adonis2(formula = HEdf_dc ~ ATTRIBUTE_blood_procedure, data = md, permutations = 999, method = "jaccard", by = NULL)
# Df SumOfSqs      R2      F Pr(>F)
# Model     1    3.762 0.04156 1.8645  0.163 #0.167 in qiime
# Residual 43   86.761 0.95844              
# Total    44   90.523 1.00000 

# adonis2(formula = HEdf_dc ~ ATTRIBUTE_blood_procedure, data = md, permutations = 999, method = "braycurtis", by = NULL)
# Df SumOfSqs      R2      F Pr(>F)   
# Model     1   0.2432 0.04866 2.1995  0.007 **
#   Residual 43   4.7547 0.95134                 
# Total    44   4.9980 1.00000  

#peripheral
md<-fread("data_files/LCMS/metadata_table/metadata_table-00000.tsv")%>%dplyr::rename(sample_name=filename)%>%
  filter(ATTRIBUTE_blood_origin=="peripheral",ATTRIBUTE_blood_procedure!="het")%>%dplyr::rename(sampleid=ATTRIBUTE_sampleID)%>%select(sampleid,everything())%>%select(-sample_name)
md$ATTRIBUTE_blood_procedure<-factor(md$ATTRIBUTE_blood_procedure,c("pre","post","bd","return"))
write.table(md,"data_files/LCMS/metadata_table/metadata_table-peripheral.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("data_files/LCMS/metadata_table/metadata_table-00000.tsv")%>%dplyr::rename(sample_name=filename)%>%
  filter(ATTRIBUTE_blood_origin=="peripheral",(ATTRIBUTE_blood_procedure=="post"|ATTRIBUTE_blood_procedure=="bd"))%>%dplyr::rename(sampleid=ATTRIBUTE_sampleID)%>%select(sampleid,everything())%>%select(-sample_name)
md$ATTRIBUTE_blood_procedure<-factor(md$ATTRIBUTE_blood_procedure,c("pre","post","bd","return"))
write.table(md,"data_files/LCMS/metadata_table/metadata_table-peripheral-postbd.txt",sep = "\t",row.names = FALSE,quote=FALSE)

perif_ord <- read_qza("diversity_metrics_LCMS/ordination-peripheral.qza")

rpca<-perif_ord$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  dplyr::rename(sampleid=SampleID)%>%
  left_join(md,by="sampleid")

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=ATTRIBUTE_blood_procedure)) +
  geom_point(alpha=1.0) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +stat_ellipse(type = "t", linetype = 2,aes(group = ATTRIBUTE_blood_procedure))+
  scale_color_manual(values=c("#52A7DC","#DF643B","#2D653D"))+
  labs(x =paste("PC1 (",round(perif_ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(perif_ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       color="TIPS visit")+ggtitle("Peripheral")+ theme(plot.title = element_text(face = "bold"))
ggsave("diversity_metrics_LCMS/SFR21_1025_PCoAvisitPerifOnlyAnnot_wellipse.tiff", plot=p,height=3, width=3.1)
ggsave("diversity_metrics_LCMS/SFR22_0512_PCoAvisitPerifOnlyAnnot_wellipse.pdf", plot=p,height=3, width=3.1)
#ggsave("diversity_metrics_LCMS/SFR21_1006_PCoAvisitPerifOnlyAnnot_braycurtis_wellipse.tiff", plot=p,height=3, width=3.1)

baseplot<-
  ggplot() +
  theme_bw() +
  labs(x =paste("PC1 (",round(perif_ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(perif_ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       fill="TIPS visit")+ggtitle("Peripheral")+
  geom_segment(data=perif_ord$data$Species %>% 
                 mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
                 top_n(5, a) %>% #keep 8 furthest away points
                 mutate(PC1=0.17*PC1, PC2=0.21*PC2)%>% # scale arrows linearly... is this ok? 
                 left_join(mtb_perif,by="FeatureID"),
               aes(x=0, xend=PC1, y=0, yend=PC2,color=feature_name),
               arrow = arrow(length = unit(0.3,"cm"))
  )

# now overlay samples
baseplot<-baseplot +
  geom_point(
    data=rpca,
    aes(x=PC1, y=PC2, fill=ATTRIBUTE_blood_procedure), 
    shape=21
  ) + scale_fill_manual(values=c("#52A7DC","#DF643B","#2D653D")) +theme(plot.title = element_text(face = "bold"))
ggsave("diversity_metrics_LCMS/SFR22_0228_PCoAvisitPerifOnly_warrows.tiff", plot=baseplot,height=4, width=19)

#10% highest/lowest ranked features -natural log
df_div<-fread("diversity_metrics_LCMS/ordination-peripheral/sample_plot_data_peripheral.tsv")%>%select(1:3)%>%
  mutate(ATTRIBUTE_blood_procedure=factor(ATTRIBUTE_blood_procedure, levels = c("pre", "post","bd")))

p<-ggplot(df_div, aes(x=ATTRIBUTE_blood_procedure, y=Current_Natural_Log_Ratio,fill=ATTRIBUTE_blood_procedure)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#52A7DC","#DF643B","#2D653D"))+scale_fill_manual(values=c("#52A7DC","#DF643B","#2D653D"))+
  labs(x="TIPS visit",y="Natural Log Ratio")+
  theme(legend.position = "none")

ggsave("diversity_metrics_LCMS/ordination-peripheral/SFR22_0228_deicode_natlog.pdf",p, height=3, width=3.5)

pairwise.wilcox.test(df_div$Current_Natural_Log_Ratio, df_div$ATTRIBUTE_blood_procedure,
                     p.adjust.method="fdr")
# pre   post 
# post 0.033 -    
# bd   0.020 2e-07

rank_feat<-fread("diversity_metrics_LCMS/ordination-peripheral/selected_features_peripheral.tsv")
#plot the lowest ranked BAs & AA
rpca<-perif_ord$data$Species %>%
  select(FeatureID, PC1, PC2) %>%
  left_join(mtb_perif,by="FeatureID")%>%
  filter(FeatureID %in% rank_feat$FeatureID)%>%
  mutate(feature_name=fct_reorder(feature_name,PC1))

p<-ggplot(data=rpca, aes(x=feature_name, y= PC1)) +
  geom_bar(stat="identity",fill="#A9A9A9")+ theme_minimal()+coord_flip()+
  labs(x="Feature",y="PC1",title="highest and lowest ranking features")
ggsave("diversity_metrics_LCMS/ordination-peripheral/SFR22_0304_topbottom10p_PC1.pdf",p, height=15, width=25)

#amino acids/bile acid ranked features -natural log
#df_div<-fread("diversity_metrics_LCMS/ordination-peripheral/sample_plot_data_bab.tsv")%>%select(1:3)%>%
df_div<-fread("diversity_metrics_LCMS/ordination-peripheral/sample_plot_data_bab_flip.tsv")%>%select(1:3)%>%
  mutate(ATTRIBUTE_blood_procedure=factor(ATTRIBUTE_blood_procedure, levels = c("pre", "post","bd")))

p<-ggplot(df_div, aes(x=ATTRIBUTE_blood_procedure, y=Current_Natural_Log_Ratio,fill=ATTRIBUTE_blood_procedure)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#52A7DC","#DF643B","#2D653D"))+scale_fill_manual(values=c("#52A7DC","#DF643B","#2D653D"))+
  labs(x="TIPS visit",y="Natural Log Ratio")+
  theme(legend.position = "none")

#ggsave("diversity_metrics_LCMS/ordination-peripheral/SFR22_0307_deicode_bab_natlog.pdf",p, height=3, width=3.5)
ggsave("diversity_metrics_LCMS/ordination-peripheral/SFR22_0718_deicode_bab_natlog.pdf",p, height=3, width=3.5)

pairwise.wilcox.test(df_div$Current_Natural_Log_Ratio, df_div$ATTRIBUTE_blood_procedure,
                     p.adjust.method="fdr")
#pre post   
#pre    post  
#post 0.0020 -     
#  bd   0.7180 0.0032

rank_feat<-fread("diversity_metrics_LCMS/ordination-peripheral/selected_features_bab.tsv")
#plot the lowest ranked BAs & Bilirubin
rpca<-perif_ord$data$Species %>%
  select(FeatureID, PC1, PC2) %>%
  left_join(mtb_perif,by="FeatureID")%>%
  filter(FeatureID %in% rank_feat$FeatureID)%>%
  mutate(feature_name=fct_reorder(feature_name,PC1))

p<-ggplot(data=rpca, aes(x=feature_name, y= PC1)) +
  geom_bar(stat="identity",fill="#A9A9A9")+ theme_minimal()+coord_flip()+#theme(axis.text.x=element_text(angle=-90))
  labs(x="Feature",y="PC1",title="Bilirubins vs. BA ranking")
ggsave("diversity_metrics_LCMS/ordination-peripheral/SFR22_0304_BilvBA_PC1.pdf",p, height=6, width=28)

HEdf_Perif<-fread("data_files/LCMS/quantification_table/quantification_table-00000.csv")%>%
  select(`row ID`, all_of(md$sampleid))
write.table(HEdf_Perif,"/data/Stephany/liverpath/data_files/LCMS/quantification_table/quantification_table-peripheral.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#HEdf_dc<-fread("diversity_metrics_LCMS/beta_deicode_lcms_peripheral/distance-matrix.tsv")%>%column_to_rownames("V1")
HEdf_dc<-fread("diversity_metrics_LCMS/beta_braycurtis_lcms_peripheral/distance-matrix.tsv")%>%column_to_rownames("V1")
HRdc_pcoa <- pcoa(HEdf_dc)
HRdc_pcoa_df<-as.data.frame(HRdc_pcoa$vectors)%>%rownames_to_column("sample_name")
md<-fread("/data/Stephany/liverpath/data_files/LCMS/metadata_table/metadata_table-peripheral.txt")%>%rename(sample_name=sampleid)%>%
  arrange(ATTRIBUTE_blood_procedure,ATTRIBUTE_pt_id)

md$ATTRIBUTE_blood_procedure[md$ATTRIBUTE_Sample=="Pre-HET"]<-"pre"
md$ATTRIBUTE_blood_procedure[grepl("Post-",md$ATTRIBUTE_Sample)]<-"post"
md$ATTRIBUTE_blood_procedure[grepl("POD-1",md$ATTRIBUTE_Sample)]<-"post"
md$ATTRIBUTE_blood_procedure<-factor(md$ATTRIBUTE_blood_procedure,c("pre","post","bd","return"))
# md$modeling<-c(rep("Test",4),rep("Train",16),rep("Test",4),rep("Train",16),rep("Test",4),rep("Train",16))
# write.table(md,"/data/Stephany/liverpath/data_files/LCMS/metadata_table/metadata_table-peripheral-TrainTest.txt",sep = "\t",row.names = FALSE,quote=FALSE)

p<-HRdc_pcoa_df %>%
  select(sample_name, Axis.1,Axis.2) %>%
  left_join(.,md, by="sample_name") %>%
  ggplot(aes(x=Axis.1, y=Axis.2, color=ATTRIBUTE_blood_procedure)) +
  geom_point(alpha=1.0) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +stat_ellipse(type = "t", linetype = 2,aes(group = ATTRIBUTE_blood_procedure))+
  scale_color_manual(values=c("#52A7DC","#DF643B","#2D653D"))+
  labs(x =paste("PC1 (",round(HRdc_pcoa$values$Relative_eig[1]*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(HRdc_pcoa$values$Relative_eig[2]*100,digits=2),"%)",sep=""),
       color="TIPS visit")+ggtitle("Peripheral")+ theme(plot.title = element_text(face = "bold"))
#ggsave("diversity_metrics_LCMS/SFR21_0929_PCoAvisitPerifOnlyAnnot_wellipse.tiff", plot=p,height=3, width=3)
ggsave("diversity_metrics_LCMS/SFR21_1006_PCoAvisitPerifOnlyAnnot_braycurtis_wellipse.tiff", plot=p,height=3, width=3.1)

md<-fread("/data/Stephany/liverpath/data_files/LCMS/metadata_table/metadata_table-peripheral.txt")%>%rename(sample_name=sampleid)%>%
  filter(ATTRIBUTE_blood_procedure=="post"|ATTRIBUTE_blood_procedure=="bd")

#HEdf_dc<-fread("diversity_metrics_LCMS/beta_deicode_lcms_peripheral/distance-matrix.tsv")%>%select("V1",contains(md$sample_name))%>%
HEdf_dc<-fread("diversity_metrics_LCMS/beta_braycurtis_lcms_peripheral/distance-matrix.tsv")%>%select("V1",contains(md$sample_name))%>%
  filter(V1 %in% md$sample_name)%>%column_to_rownames("V1")
adonis2(HEdf_dc ~ ATTRIBUTE_blood_procedure, data = md, permutations = 999, method="braycurtis",by = NULL)

#pre vs post
# adonis2(formula = HEdf_dc ~ ATTRIBUTE_blood_procedure, data = md, permutations = 999, method = "jaccard", by = NULL)
# Df SumOfSqs      R2      F Pr(>F)  
# Model     1    6.998 0.10208 4.3199  0.017 * #0.013 qiime   q=0.0195
#   Residual 38   61.553 0.89792                
# Total    39   68.551 1.00000 

# adonis2(formula = HEdf_dc ~ ATTRIBUTE_blood_procedure, data = md, permutations = 999, method = "braycurtis", by = NULL)
# Df SumOfSqs      R2      F Pr(>F)   
# Model     1   0.2785 0.06586 2.6791  0.002 **
#   Residual 38   3.9507 0.93414                 
# Total    39   4.2292 1.00000 

#pre vs. bd
# adonis2(formula = HEdf_dc ~ ATTRIBUTE_blood_procedure, data = md, permutations = 999, method = "jaccard", by = NULL)
# Df SumOfSqs      R2      F Pr(>F)
# Model     1    4.323 0.05082 2.0347  0.117 #0.141 q is the same
# Residual 38   80.728 0.94918              
# Total    39   85.050 1.00000 

# adonis2(formula = HEdf_dc ~ ATTRIBUTE_blood_procedure, data = md, permutations = 999, method = "braycurtis", by = NULL)
# Df SumOfSqs      R2      F Pr(>F)  
# Model     1   0.2108 0.04318 1.7147  0.054 . 
# Residual 38   4.6720 0.95682                
# Total    39   4.8828 1.00000 

#post vs bd
# adonis2(formula = HEdf_dc ~ ATTRIBUTE_blood_procedure, data = md, permutations = 999, method = "jaccard", by = NULL)
# Df SumOfSqs      R2      F Pr(>F)    
# Model     1   22.197 0.28896 15.443  0.001 *** #0.001 q=0.0030
#   Residual 38   54.619 0.71104                  
# Total    39   76.816 1.00000  

# adonis2(formula = HEdf_dc ~ ATTRIBUTE_blood_procedure, data = md, permutations = 999, method = "braycurtis", by = NULL)
# Df SumOfSqs      R2      F Pr(>F)    
# Model     1   0.5197 0.11856 5.1111  0.001 ***
#   Residual 38   3.8639 0.88144                  
# Total    39   4.3836 1.00000

#peripheral post - BD
md<-fread("data_files/LCMS/metadata_table/metadata_table-00000.tsv")%>%dplyr::rename(sample_name=filename)%>%
  filter(ATTRIBUTE_blood_origin=="peripheral",(ATTRIBUTE_blood_procedure=="post"|ATTRIBUTE_blood_procedure=="bd"))%>%dplyr::rename(sampleid=ATTRIBUTE_sampleID)%>%select(sampleid,everything())%>%select(-sample_name)
md$ATTRIBUTE_blood_procedure<-factor(md$ATTRIBUTE_blood_procedure,c("pre","post","bd","return"))
write.table(md,"data_files/LCMS/metadata_table/metadata_table-peripheral-postbd.txt",sep = "\t",row.names = FALSE,quote=FALSE)

HEdf_Perif<-fread("data_files/LCMS/quantification_table/quantification_table-00000.csv")%>%
  select(`row ID`, all_of(md$sampleid))
write.table(HEdf_Perif,"data_files/LCMS/quantification_table/quantification_table-peripheral-postbd.txt",sep = "\t",row.names = FALSE,quote=FALSE)

perif_ord <- read_qza("diversity_metrics_LCMS/ordination-peripheral-postbd.qza")

rpca<-perif_ord$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  dplyr::rename(sampleid=SampleID)%>%
  left_join(md,by="sampleid")

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=ATTRIBUTE_blood_procedure)) +
  geom_point(alpha=1.0) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +stat_ellipse(type = "t", linetype = 2,aes(group = ATTRIBUTE_blood_procedure))+
  scale_color_manual(values=c("#DF643B","#2D653D"))+
  labs(x =paste("PC1 (",round(perif_ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(perif_ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       color="TIPS visit")+ggtitle("Peripheral")+ theme(plot.title = element_text(face = "bold"))
ggsave("diversity_metrics_LCMS/SFR22_0718_PCoAvisitPerifPostBDAnnot_wellipse.pdf", plot=p,height=3, width=3.1)

#just post
md<-fread("data_files/LCMS/metadata_table/metadata_table-00000.tsv")

# md$ATTRIBUTE_blood_procedure[md$ATTRIBUTE_Sample=="Pre-HET"]<-"pre"
# md$ATTRIBUTE_blood_procedure[grepl("Post-",md$ATTRIBUTE_Sample)]<-"post"
# md$ATTRIBUTE_blood_procedure[grepl("POD-1",md$ATTRIBUTE_Sample)]<-"post"
md$ATTRIBUTE_blood_procedure<-factor(md$ATTRIBUTE_blood_procedure,c("pre","post","bd","return"))
md$ATTRIBUTE_blood_origin<-factor(md$ATTRIBUTE_blood_origin,c("peripheral","hepatic"))
md<-md%>%dplyr::rename(sample_name=filename)%>%
  filter(ATTRIBUTE_blood_procedure=="post")%>%dplyr::rename(sampleid=ATTRIBUTE_sampleID)%>%select(sampleid,everything())%>%select(-sample_name)%>%
  mutate(WPostTIPS_HE=ifelse(ATTRIBUTE_Worst_PostTIPS_HE==3|ATTRIBUTE_Worst_PostTIPS_HE==2,"2+",ATTRIBUTE_Worst_PostTIPS_HE))
write.table(md,"/data/Stephany/liverpath/data_files/LCMS/metadata_table/metadata_table-post.txt",sep = "\t",row.names = FALSE,quote=FALSE)

grd_ord <- read_qza("diversity_metrics_LCMS/ordination-post.qza")

rpca<-grd_ord$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  dplyr::rename(sampleid=SampleID)%>%
  left_join(md,by="sampleid")

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=WPostTIPS_HE,shape=ATTRIBUTE_blood_origin)) +
  geom_point(alpha=1.0) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +stat_ellipse(type = "t", linetype = 2,aes(group = WPostTIPS_HE))+
  scale_shape_manual(values=c(15,16,17),name="Blood Origin") +
  scale_color_manual(values=c("#39B54A", "#283891", "#EF3E36"))+
  labs(x =paste("PC1 (",round(grd_ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(grd_ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       color="HE grade")+ggtitle("PostTIPS") + theme(plot.title = element_text(face = "bold"), legend.position = "left")
ggsave("diversity_metrics_LCMS/SFR21_1025_PCoAvisitPostOnlyHEgrade_wellipse.tiff", plot=p,height=3, width=4.5)
ggsave("diversity_metrics_LCMS/SFR22_0512_PCoAvisitPostOnlyHEgrade_wellipse.pdf", plot=p,height=3, width=4.5)

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=ATTRIBUTE_blood_origin)) +
  geom_point(alpha=1.0) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +stat_ellipse(type = "t", linetype = 2,aes(group = ATTRIBUTE_blood_origin))+
  scale_shape_manual(values=c(15,16,17),name="Blood Origin") +
  scale_color_manual(values=c("#2E3192","#EE207C")) +
  labs(x =paste("PC1 (",round(grd_ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(grd_ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       color="Blood Origin")+ggtitle("PostTIPS") + theme(plot.title = element_text(face = "bold"), legend.position = "top")
ggsave("diversity_metrics_LCMS/SFR22_0426_PCoAvisitPostOnlyBldOrg_wellipse.tiff", plot=p,height=3, width=4.5)
#ggsave("diversity_metrics_LCMS/SFR21_1006_PCoAvisitPostOnlyHEgrade_braycurtis_wellipse.tiff", plot=p,height=3, width=4.5)

baseplot<-
  ggplot() +
  theme_bw() +
  labs(x =paste("PC1 (",round(grd_ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(grd_ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       fill="TIPS visit")+ggtitle("Peripheral")+
  geom_segment(data=grd_ord$data$Species %>% 
                 mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
                 top_n(5, a) %>% #keep 8 furthest away points
                 mutate(PC1=0.26*PC1, PC2=0.18*PC2)%>% # scale arrows linearly... is this ok? 
                 left_join(mtb_perif,by="FeatureID"),
               aes(x=0, xend=PC1, y=0, yend=PC2,color=feature_name),
               arrow = arrow(length = unit(0.2,"cm"))
  )

# now overlay samples
baseplot<-baseplot +
  geom_point(
    data=rpca,
    aes(x=PC1, y=PC2, fill=WPostTIPS_HE), 
    shape=21
  ) + scale_fill_manual(values=c("#DC4D40","#754B27","#60B358")) +theme(plot.title = element_text(face = "bold"))
ggsave("diversity_metrics_LCMS/SFR22_0228_PCoAvisitPostOnlyHEgrade_warrows.tiff", plot=baseplot,height=4, width=20)

#10% highest/lowest ranked features -natural log
##based on HE grade
df_div<-fread("diversity_metrics_LCMS/ordination-post/sample_plot_data_post.tsv")%>%
  mutate(ATTRIBUTE_Worst_PostTIPS_HE_mod=factor(ATTRIBUTE_Worst_PostTIPS_HE_mod, levels = c("0", "1","2+")),
         ATTRIBUTE_blood_origin=factor(ATTRIBUTE_blood_origin,levels=c("peripheral","hepatic")))

p<-ggplot(df_div, aes(x=ATTRIBUTE_Worst_PostTIPS_HE_mod, y=Current_Natural_Log_Ratio,fill=ATTRIBUTE_Worst_PostTIPS_HE_mod)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) + facet_wrap(~ATTRIBUTE_blood_origin,ncol=2,scales="free_y") +
  theme_minimal()+ scale_color_manual(values=c("#DC4D40","#754B27","#60B358"))+scale_fill_manual(values=c("#DC4D40","#754B27","#60B358"))+
  labs(x="HE grade",y="Natural Log Ratio")+
  theme(legend.position = "none")

ggsave("diversity_metrics_LCMS/ordination-post/SFR22_0301_SFR22_0228_deicode_natlog.pdf",p, height=3, width=5)

per<-df_div%>%filter(ATTRIBUTE_blood_origin=="hepatic")
pairwise.wilcox.test(df_div$Current_Natural_Log_Ratio, df_div$ATTRIBUTE_Worst_PostTIPS_HE_mod,
                     p.adjust.method="fdr")
#per
#0 1
# 1  1 -
#   2+ 1 1

#hep
# 0    1   
# 1  0.87 -   
#   2+ 0.87 0.87

#comb
# 0    1   
# 1  0.76 -   
#   2+ 0.76 0.76

#based on per vs hep
df_div<-fread("diversity_metrics_LCMS/ordination-post/sample_plot_data_post.tsv")%>%
  mutate(ATTRIBUTE_Worst_PostTIPS_HE_mod=factor(ATTRIBUTE_Worst_PostTIPS_HE_mod, levels = c("0", "1","2+")),
         ATTRIBUTE_blood_origin=factor(ATTRIBUTE_blood_origin,levels=c("peripheral","hepatic")))

p<-ggplot(df_div, aes(x=ATTRIBUTE_blood_origin, y=Current_Natural_Log_Ratio,fill=ATTRIBUTE_blood_origin)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#2E3192","#EE207C"))+scale_fill_manual(values=c("#2E3192","#EE207C"))+
  labs(x="Blood Origin",y="Natural Log Ratio")+
  theme(legend.position = "none")

ggsave("diversity_metrics_LCMS/ordination-post/SFR22_0301_deicode_bldorg_natlog.pdf",p, height=3, width=2.5)

pairwise.wilcox.test(df_div$Current_Natural_Log_Ratio, df_div$ATTRIBUTE_blood_origin,
                     p.adjust.method="fdr")

# peripheral
# hepatic 0.72

HEdf_Post<-fread("data_files/LCMS/quantification_table/quantification_table-00000.csv")%>%
  select(`row ID`, all_of(md$sampleid))
write.table(HEdf_Post,"/data/Stephany/liverpath/data_files/LCMS/quantification_table/quantification_table-post.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("/data/Stephany/liverpath/data_files/LCMS/metadata_table/metadata_table-post.txt")%>%arrange(WPostTIPS_HE)
md$modeling<-c(rep("Test",2),rep("Train",6),rep("Test",2),rep("Train",17),rep("Test",2),rep("Train",14))
write.table(md,"/data/Stephany/liverpath/data_files/LCMS/metadata_table/metadata_table-post-TrainTest.txt",sep = "\t",row.names = FALSE,quote=FALSE)

HEdf_dc<-read_qza("diversity_metrics_LCMS/ctf-results/state_biplot.qza")


#HEdf_dc<-fread("diversity_metrics_LCMS/beta_deicode_lcms_post/distance-matrix.tsv")%>%column_to_rownames("V1")
#HEdf_dc<-fread("diversity_metrics_LCMS/beta_braycurtis_lcms_post/distance-matrix.tsv")%>%column_to_rownames("V1")
HRdc_pcoa <- pcoa(HEdf_dc)
HRdc_pcoa_df<-as.data.frame(HRdc_pcoa$vectors)%>%rownames_to_column("sampleid")

p<-HRdc_pcoa_df %>%
  select(sampleid, Axis.1,Axis.2) %>%
  left_join(.,md, by="sampleid") %>%
  ggplot(aes(x=Axis.1, y=Axis.2, color=WPostTIPS_HE,shape=ATTRIBUTE_blood_origin)) +
  geom_point(alpha=1.0) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +stat_ellipse(type = "t", linetype = 2,aes(group = WPostTIPS_HE))+
  scale_shape_manual(values=c(15,16,17),name="Blood Origin") +
  scale_color_manual(values=c("#DC4D40","#754B27","#60B358"))+
  labs(x =paste("PC1 (",round(HRdc_pcoa$values$Relative_eig[1]*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(HRdc_pcoa$values$Relative_eig[2]*100,digits=2),"%)",sep=""),
       color="HE grade")+ggtitle("PostTIPS") + theme(plot.title = element_text(face = "bold"), legend.position = "left")
#ggsave("diversity_metrics_LCMS/SFR21_0930_PCoAvisitPostOnlyHEgrade_wellipse.tiff", plot=p,height=3, width=4.5)
ggsave("diversity_metrics_LCMS/SFR21_1006_PCoAvisitPostOnlyHEgrade_braycurtis_wellipse.tiff", plot=p,height=3, width=4.5)

md<-fread("/data/Stephany/liverpath/data_files/LCMS/metadata_table/metadata_table-post.txt")%>%rename(sample_name=sampleid)%>%
  filter(WPostTIPS_HE=="1"|WPostTIPS_HE=="2+")

#HEdf_dc<-fread("diversity_metrics_LCMS/beta_deicode_lcms_post/distance-matrix.tsv")%>%select("V1",contains(md$sample_name))%>%
HEdf_dc<-fread("diversity_metrics_LCMS/beta_braycurtis_lcms_post/distance-matrix.tsv")%>%select("V1",contains(md$sample_name))%>%
  filter(V1 %in% md$sample_name)%>%column_to_rownames("V1")
adonis2(HEdf_dc ~ WPostTIPS_HE, data = md, permutations = 999, method="braycurtis",by = NULL)

#0 vs 2+
# adonis2(formula = HEdf_dc ~ WPostTIPS_HE, data = md, permutations = 999, method = "jaccard", by = NULL)
# Df SumOfSqs      R2      F Pr(>F)  
# Model     1    6.167 0.11377 2.8244  0.069 . #0.071 q=0.108
# Residual 22   48.037 0.88623                
# Total    23   54.204 1.00000 

# adonis2(formula = HEdf_dc ~ WPostTIPS_HE, data = md, permutations = 999, method = "braycurtis", by = NULL)
# Df SumOfSqs      R2      F Pr(>F)
# Model     1  0.10053 0.04181 0.9601  0.476
# Residual 22  2.30376 0.95819              
# Total    23  2.40429 1.00000

#1 vs. 2+
# adonis2(formula = HEdf_dc ~ WPostTIPS_HE, data = md, permutations = 999, method = "jaccard", by = NULL)
# Df SumOfSqs      R2      F Pr(>F)  
# Model     1    6.288 0.07784 2.7856  0.058 . #0.072 q=0.108
# Residual 33   74.486 0.92216                
# Total    34   80.774 1.00000   

# adonis2(formula = HEdf_dc ~ WPostTIPS_HE, data = md, permutations = 999, method = "braycurtis", by = NULL)
# Df SumOfSqs      R2      F Pr(>F)  
# Model     1  0.13574 0.04307 1.4852   0.06 .
# Residual 33  3.01597 0.95693                
# Total    34  3.15171 1.00000 

#0 vs. 1
# adonis2(formula = HEdf_dc ~ WPostTIPS_HE, data = md, permutations = 999, method = "jaccard", by = NULL)
# Df SumOfSqs      R2      F Pr(>F)
# Model     1    0.309 0.00923 0.2329  0.842 #0.820
# Residual 25   33.130 0.99077              
# Total    26   33.438 1.00000 

# adonis2(formula = HEdf_dc ~ WPostTIPS_HE, data = md, permutations = 999, method = "braycurtis", by = NULL)
# Df SumOfSqs      R2      F Pr(>F)
# Model     1  0.07144 0.03611 0.9366  0.509
# Residual 25  1.90680 0.96389              
# Total    26  1.97824 1.00000 

#just pre
md<-fread("data_files/LCMS/metadata_table/metadata_table-00000.tsv")

# md$ATTRIBUTE_blood_procedure[md$ATTRIBUTE_Sample=="Pre-HET"]<-"pre"
# md$ATTRIBUTE_blood_procedure[grepl("Post-",md$ATTRIBUTE_Sample)]<-"post"
# md$ATTRIBUTE_blood_procedure[grepl("POD-1",md$ATTRIBUTE_Sample)]<-"post"
md$ATTRIBUTE_blood_procedure<-factor(md$ATTRIBUTE_blood_procedure,c("pre","post","bd","return"))
md$ATTRIBUTE_blood_origin<-factor(md$ATTRIBUTE_blood_origin,c("peripheral","hepatic"))
md<-md%>%dplyr::rename(sample_name=filename)%>%
  filter(ATTRIBUTE_blood_procedure=="pre")%>%dplyr::rename(sampleid=ATTRIBUTE_sampleID)%>%select(sampleid,everything())%>%select(-sample_name)%>%
  mutate(WPostTIPS_HE=ifelse(ATTRIBUTE_Worst_PostTIPS_HE==3|ATTRIBUTE_Worst_PostTIPS_HE==2,"2+",ATTRIBUTE_Worst_PostTIPS_HE))
write.table(md,"data_files/LCMS/metadata_table/metadata_table-pre.txt",sep = "\t",row.names = FALSE,quote=FALSE)

HEdf_Pre<-fread("data_files/LCMS/quantification_table/quantification_table-00000.csv")%>%
  select(`row ID`, all_of(md$sampleid))
write.table(HEdf_Pre,"data_files/LCMS/quantification_table/quantification_table-pre.txt",sep = "\t",row.names = FALSE,quote=FALSE)

grd_ord <- read_qza("diversity_metrics_LCMS/ordination-pre.qza")

rpca<-grd_ord$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  dplyr::rename(sampleid=SampleID)%>%
  left_join(md,by="sampleid")

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=factor(ATTRIBUTE_Worst_PreTIPS_HE),shape=ATTRIBUTE_blood_origin)) +
  geom_point(alpha=1.0) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +stat_ellipse(type = "t", linetype = 2,aes(group = factor(ATTRIBUTE_Worst_PreTIPS_HE)))+
  scale_shape_manual(values=c(15,16,17),name="Blood Origin") +
  scale_color_manual(values=c("#DC4D40","#754B27","#60B358"))+
  labs(x =paste("PC1 (",round(grd_ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(grd_ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       color="HE grade")+ggtitle("PreTIPS") + theme(plot.title = element_text(face = "bold"), legend.position = "left")
ggsave("diversity_metrics_LCMS/SFR22_0304_PCoAvisitPreOnlyHEgrade_wellipse.tiff", plot=p,height=3, width=4.5)


p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=ATTRIBUTE_blood_origin)) +
  geom_point(alpha=1.0) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +stat_ellipse(type = "t", linetype = 2,aes(group = ATTRIBUTE_blood_origin))+
  scale_shape_manual(values=c(15,16,17),name="Blood Origin") +
  scale_color_manual(values=c("#2E3192","#EE207C")) +
  labs(x =paste("PC1 (",round(grd_ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(grd_ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       color="Blood Origin")+ggtitle("PreTIPS") + theme(plot.title = element_text(face = "bold"), legend.position = "top")
ggsave("diversity_metrics_LCMS/SFR22_0304_PCoAvisitPreOnlyBldOrg_wellipse.tiff", plot=p,height=3.5, width=4)
#ggsave("diversity_metrics_LCMS/SFR21_1006_PCoAvisitPostOnlyHEgrade_braycurtis_wellipse.tiff", plot=p,height=3, width=4.5)

baseplot<-
  ggplot() +
  theme_bw() +
  labs(x =paste("PC1 (",round(grd_ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(grd_ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       fill="TIPS visit")+ggtitle("Peripheral")+
  geom_segment(data=grd_ord$data$Species %>% 
                 mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
                 top_n(5, a) %>% #keep 8 furthest away points
                 mutate(PC1=0.26*PC1, PC2=0.18*PC2)%>% # scale arrows linearly... is this ok? 
                 left_join(mtb_perif,by="FeatureID"),
               aes(x=0, xend=PC1, y=0, yend=PC2,color=feature_name),
               arrow = arrow(length = unit(0.2,"cm"))
  )

# now overlay samples
baseplot<-baseplot +
  geom_point(
    data=rpca,
    aes(x=PC1, y=PC2, fill=WPostTIPS_HE), 
    shape=21
  ) + scale_fill_manual(values=c("#DC4D40","#754B27","#60B358")) +theme(plot.title = element_text(face = "bold"))
ggsave("diversity_metrics_LCMS/SFR22_0304_PCoAvisitPreOnlyBldOrg_warrows.tiff", plot=baseplot,height=4, width=20)

#10% highest/lowest ranked features -natural log
##based on HE grade
df_div<-fread("diversity_metrics_LCMS/ordination-pre/sample_plot_data_pre.tsv")%>%
  mutate(ATTRIBUTE_Worst_PreTIPS_HE_mod=factor(ATTRIBUTE_Worst_PreTIPS_HE_mod, levels = c("0", "1","2+")),
         ATTRIBUTE_blood_origin=factor(ATTRIBUTE_blood_origin,levels=c("peripheral","hepatic")))

p<-ggplot(df_div, aes(x=ATTRIBUTE_Worst_PreTIPS_HE_mod, y=Current_Natural_Log_Ratio,fill=ATTRIBUTE_Worst_PreTIPS_HE_mod)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) + facet_wrap(~ATTRIBUTE_blood_origin,ncol=2,scales="free_y") +
  theme_minimal()+ scale_color_manual(values=c("#DC4D40","#754B27","#60B358"))+scale_fill_manual(values=c("#DC4D40","#754B27","#60B358"))+
  labs(x="HE grade",y="Natural Log Ratio")+
  theme(legend.position = "none")

ggsave("diversity_metrics_LCMS/ordination-pre/SFR22_0304_deicode_natlog.pdf",p, height=3, width=4)

per<-df_div%>%filter(ATTRIBUTE_blood_origin=="hepatic")
pairwise.wilcox.test(df_div$Current_Natural_Log_Ratio, df_div$ATTRIBUTE_Worst_PreTIPS_HE_mod,
                     p.adjust.method="fdr")
#per
#0    1   
#1  0.55 -   
#  2+ 0.55 0.60

#hep
#0    1   
#1  0.78 -   
#  2+ 0.55 0.55

#comb
#0    1   
#1  0.43 -   
#  2+ 0.23 0.25

#based on per vs hep
df_div<-fread("diversity_metrics_LCMS/ordination-pre/sample_plot_data_pre.tsv")%>%
  mutate(ATTRIBUTE_Worst_PreTIPS_HE_mod=factor(ATTRIBUTE_Worst_PreTIPS_HE_mod, levels = c("0", "1")),
         ATTRIBUTE_blood_origin=factor(ATTRIBUTE_blood_origin,levels=c("peripheral","hepatic")))%>%
  filter(!is.na(ATTRIBUTE_Worst_PreTIPS_HE_mod))

p<-ggplot(df_div, aes(x=ATTRIBUTE_blood_origin, y=Current_Natural_Log_Ratio,fill=ATTRIBUTE_blood_origin)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#2E3192","#EE207C"))+scale_fill_manual(values=c("#2E3192","#EE207C"))+
  labs(x="Blood Origin",y="Natural Log Ratio")+
  theme(legend.position = "none")

ggsave("diversity_metrics_LCMS/ordination-pre/SFR22_0304_deicode_bldorg_natlog.pdf",p, height=3, width=2.5)

pairwise.wilcox.test(df_div$Current_Natural_Log_Ratio, df_div$ATTRIBUTE_blood_origin,
                     p.adjust.method="fdr")

#peripheral
#hepatic 0.046

#hepatic
md<-fread("/mnt/zarrinpar/scratch/sfloresr/HE/HE_metab/liverpath/data_files/LCMS/metadata_table/metadata_table-00000.tsv")%>%rename(sample_name=filename)%>%
  filter(ATTRIBUTE_blood_origin=="hepatic")%>%rename(sampleid=ATTRIBUTE_sampleID)%>%select(sampleid,everything())%>%select(-sample_name)
md$ATTRIBUTE_blood_procedure<-factor(md$ATTRIBUTE_blood_procedure,c("pre","post","bd","return"))
write.table(md,"data_files/LCMS/metadata_table/metadata_table-hepatic.txt",sep = "\t",row.names = FALSE,quote=FALSE)

HEdf_hepatic<-fread("data_files/LCMS/quantification_table/quantification_table-00000.csv")%>%
  select(`row ID`, all_of(md$sampleid))
write.table(HEdf_hepatic,"/data/Stephany/liverpath/data_files/LCMS/quantification_table/quantification_table-hepatic.txt",sep = "\t",row.names = FALSE,quote=FALSE)


HEdf_hep<-fread("data_files/LCMS/quantification_table/quantification_table-00000.csv")%>%
  select(`row ID`, all_of(md$sampleid))
write.table(HEdf_hep,"data_files/LCMS/quantification_table/quantification_table-hepatic-rmoutlier.txt",sep = "\t",row.names = FALSE,quote=FALSE)

hep_ord <- read_qza("diversity_metrics_LCMS/ordination-hepatic.qza")

rpca<-hep_ord$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  rename(sampleid=SampleID)%>%
  left_join(md,by="sampleid")

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=ATTRIBUTE_blood_procedure)) +
  geom_point(alpha=1.0, shape=15) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +stat_ellipse(type = "t", linetype = 2,aes(group = ATTRIBUTE_blood_procedure))+
  scale_color_manual(values=c("#52A7DC","#DF643B","#2D653D"))+
  labs(x =paste("PC1 (",round(hep_ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(hep_ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       color="TIPS visit")+ggtitle("Hepatic")+ theme(plot.title = element_text(face = "bold"))
ggsave("diversity_metrics_LCMS/SFR21_1025_PCoAvisitHVOnlyAnnot_wellipse.tiff", plot=p,height=3, width=3)
ggsave("diversity_metrics_LCMS/SFR22_0512_PCoAvisitHVOnlyAnnot_wellipse.pdf", plot=p,height=3, width=3)
#ggsave("diversity_metrics_LCMS/SFR21_1006_PCoAvisitHVOnlyAnnot_braycurtis_wellipse.tiff", plot=p,height=3, width=3)


#10% highest/lowest ranked features -natural log
df_div<-fread("diversity_metrics_LCMS/ordination-hepatic/sample_plot_data_hepatic.tsv")%>%select(1:3)%>%
  mutate(ATTRIBUTE_blood_procedure=factor(ATTRIBUTE_blood_procedure, levels = c("pre", "post","bd")))

p<-ggplot(df_div, aes(x=ATTRIBUTE_blood_procedure, y=Current_Natural_Log_Ratio,fill=ATTRIBUTE_blood_procedure)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#52A7DC","#DF643B","#2D653D"))+scale_fill_manual(values=c("#52A7DC","#DF643B","#2D653D"))+
  scale_shape_manual(values = c(15)) +
  labs(x="TIPS visit",y="Natural Log Ratio")+
  theme(legend.position = "none")

ggsave("diversity_metrics_LCMS/ordination-hepatic/SFR22_0301_deicode_natlog.pdf",p, height=3, width=2.3)

pairwise.wilcox.test(df_div$Current_Natural_Log_Ratio, df_div$ATTRIBUTE_blood_procedure,
                     p.adjust.method="fdr")
#   pre 
# post 0.82


#HEdf_dc<-fread("diversity_metrics_LCMS/beta_deicode_lcms_hepatic/distance-matrix.tsv")%>%column_to_rownames("V1")
HEdf_dc<-fread("diversity_metrics_LCMS/beta_braycurtis_lcms_hepatic/distance-matrix.tsv")%>%column_to_rownames("V1")
HRdc_pcoa <- pcoa(HEdf_dc)
HRdc_pcoa_df<-as.data.frame(HRdc_pcoa$vectors)%>%rownames_to_column("sample_name")
md<-fread("/data/Stephany/liverpath/data_files/LCMS/metadata_table/metadata_table-hepatic.txt")%>%rename(sample_name=sampleid)

md$ATTRIBUTE_blood_procedure[md$ATTRIBUTE_Sample=="Pre-HET"]<-"pre"
md$ATTRIBUTE_blood_procedure[grepl("Post-",md$ATTRIBUTE_Sample)]<-"post"
md$ATTRIBUTE_blood_procedure[grepl("POD-1",md$ATTRIBUTE_Sample)]<-"post"
md$ATTRIBUTE_blood_procedure<-factor(md$ATTRIBUTE_blood_procedure,c("pre","post","bd","return"))

p<-HRdc_pcoa_df %>%
  select(sample_name, Axis.1,Axis.2) %>%
  left_join(.,md, by="sample_name") %>%
  ggplot(aes(x=Axis.1, y=Axis.2, color=ATTRIBUTE_blood_procedure)) +
  geom_point(alpha=1.0, shape=15) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +stat_ellipse(type = "t", linetype = 2,aes(group = ATTRIBUTE_blood_procedure))+
  scale_color_manual(values=c("#52A7DC","#DF643B","#2D653D"))+
  labs(x =paste("PC1 (",round(HRdc_pcoa$values$Relative_eig[1]*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(HRdc_pcoa$values$Relative_eig[2]*100,digits=2),"%)",sep=""),
       color="TIPS visit")+ggtitle("Hepatic")+ theme(plot.title = element_text(face = "bold"))
#ggsave("diversity_metrics_LCMS/SFR21_0824_PCoAvisitHVOnlyAnnot_wellipse.tiff", plot=p,height=3, width=3)
ggsave("diversity_metrics_LCMS/SFR21_1006_PCoAvisitHVOnlyAnnot_braycurtis_wellipse.tiff", plot=p,height=3, width=3)

adonis2(HEdf_dc ~ ATTRIBUTE_blood_procedure, data = md, permutations = 999, method="braycurtis",by = NULL)

# adonis2(formula = HEdf_dc ~ ATTRIBUTE_blood_procedure, data = md, permutations = 999, method = "jaccard", by = NULL)
# Df SumOfSqs      R2      F Pr(>F)
# Model     1    0.273 0.00339 0.1293  0.877 #0.878
# Residual 38   80.284 0.99661              
# Total    39   80.557 1.00000 

# adonis2(formula = HEdf_dc ~ ATTRIBUTE_blood_procedure, data = md, permutations = 999, method = "braycurtis", by = NULL)
# Df SumOfSqs      R2      F Pr(>F)
# Model     1   0.1148 0.03013 1.1804  0.241 
# Residual 38   3.6970 0.96987              
# Total    39   3.8119 1.00000

#hepatic removing plasma sample (1b)
md<-fread("/mnt/zarrinpar/scratch/sfloresr/HE/HE_metab/liverpath/data_files/LCMS/metadata_table/metadata_table-00000.tsv")%>%rename(sample_name=filename)%>%
  filter(ATTRIBUTE_blood_origin=="hepatic" & ATTRIBUTE_pt_id!="1_b")%>%rename(sampleid=ATTRIBUTE_sampleID)%>%select(sampleid,everything())%>%select(-sample_name)
md$ATTRIBUTE_blood_procedure<-factor(md$ATTRIBUTE_blood_procedure,c("pre","post","bd","return"))
write.table(md,"data_files/LCMS/metadata_table/metadata_table-hepatic-rmoutlier.txt",sep = "\t",row.names = FALSE,quote=FALSE)

HEdf_hep<-fread("data_files/LCMS/quantification_table/quantification_table-00000.csv")%>%
  select(`row ID`, all_of(md$sampleid))
write.table(HEdf_hep,"data_files/LCMS/quantification_table/quantification_table-hepatic-rmoutlier.txt",sep = "\t",row.names = FALSE,quote=FALSE)

hep_ord <- read_qza("diversity_metrics_LCMS/ordination-hepatic-rmoutlier.qza")

rpca<-hep_ord$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  rename(sampleid=SampleID)%>%
  left_join(md,by="sampleid")

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=ATTRIBUTE_blood_procedure)) +
  geom_point(alpha=1.0, shape=15) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +stat_ellipse(type = "t", linetype = 2,aes(group = ATTRIBUTE_blood_procedure))+
  scale_color_manual(values=c("#52A7DC","#DF643B","#2D653D"))+
  labs(x =paste("PC1 (",round(hep_ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(hep_ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       color="TIPS visit")+ggtitle("Hepatic")+ theme(plot.title = element_text(face = "bold"))
ggsave("diversity_metrics_LCMS/SFR22_0718_PCoAvisitHVOnlyrmoutlierAnnot_wellipse.pdf", plot=p,height=3, width=3)

##################################################

md<-fread("/data/Stephany/liverpath/data_files/LCMS/metadata_table/metadata_table-00000.tsv")%>%rename(sample_name=filename)%>%
  select(sample_name,everything())
# md$ATTRIBUTE_blood_procedure[md$ATTRIBUTE_Sample=="Pre-HET"]<-"pre"
# md$ATTRIBUTE_blood_procedure[grepl("Post-",md$ATTRIBUTE_Sample)]<-"post"
# md$ATTRIBUTE_blood_procedure[grepl("POD-1",md$ATTRIBUTE_Sample)]<-"post"
md$ATTRIBUTE_blood_procedure<-factor(md$ATTRIBUTE_blood_procedure,c("pre","post","bd","return"))

write.table(md,"/data/Stephany/liverpath/data_files/LCMS/metadata_table/metadata_table-00000-clean.txt",sep = "\t",row.names = FALSE,quote=FALSE)

ord <- read_qza("diversity_metrics_LCMS/ordination.qza")

rpca<-ord$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, shape=ATTRIBUTE_blood_origin, color=ATTRIBUTE_blood_procedure)) +
  geom_point(alpha=1.0) + #alpha controls transparency and helps when points are overlapping
  theme_classic() +stat_ellipse(type = "t", linetype = 2,aes(group = ATTRIBUTE_blood_procedure))+
  scale_shape_manual(values=c(15,16,17),name="Blood Origin") +
  scale_color_manual(values=c("#52A7DC","#DF643B","#2D653D"))+
  labs(x =paste("PC1 (",round(ord$data$ProportionExplained[1]*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained[2]*100,digits=2),"%)",sep=""),
       color="TIPS visit")
ggsave("diversity_metrics_LCMS/SFR21_1025_PCoAvisitBldOrgAnnot_wellipse.tiff", plot=p,height=4, width=5.5)

##################################################

mtb_divers<-fread("diff_HE_analysis/SFR21_1025_HEmtbdata_allwmd.txt")%>%select(rowID,Compound_Name)%>%unique()%>%
  mutate(annotated=ifelse(is.na(Compound_Name),"Unannotated","Annotated"))%>%group_by(annotated)%>%summarise(n=n())%>%
  mutate(mtbs="mtbs",pct=n/595)
  
#create a barplot summarizing the composition--dont run
p<-ggplot(mtb_divers, aes(fill=fct_rev(annotated), y=n, x=mtbs)) + geom_bar(position="stack", stat="identity") +
  geom_text(aes(y =c(68.5,350), label = paste(n,"\n(",signif(pct*100,digits=4),"%)")), vjust = 0.7, colour = "black") + 
  scale_fill_manual(values=c('grey70','grey34'))+
  theme_void()+coord_flip()+guides(fill=guide_legend(title="n=595 mtbs"))

ggsave("diff_HE_analysis/SFR21_1025_lcmsHE_compbarplot_nounclass.tiff", plot=p,height=0.9, width=4.5)
##################################################
setwd("/data/Stephany/liverpath/")
spectral_counts<-(read_qza("qtree_analysis/feature-table-hashed.qza")$data)%>%as.data.frame()%>%
  select(-`row_m/z`,-rt_time)%>%as.matrix()
spectral_counts[spectral_counts==0] <- 1
feature_table<-fread("/data/Stephany/liverpath/qtree_analysis/classified-merged-feature-data/feature_data.tsv")%>%column_to_rownames("V1")%>%
  select("kingdom","superclass","class","subclass","direct_parent","smiles")%>%as.matrix()

#setwd("/mnt/genomics/Stephany/liverpath/LCMS/")
# md_pre<-fread("metadata_table/SFR21_0513_metadata_table_updatedTIPS_long_v2.txt")%>%
#   select(sampleID,blood_draw,ATTRIBUTE_blood_procedure,TIPS_HE,ammonia,ATTRIBUTE_diagnosis_short)%>%
#   replace_na(list(blood_draw = "none", ATTRIBUTE_blood_procedure = "none",TIPS_HE="none",ammonia="0",ATTRIBUTE_diagnosis_short="none"))%>%
#   as.matrix()
#write.table(md_pre,"metadata_table/SFR21_0514_metadata_table_updatedTIPS_long_v3.txt",sep = "\t",row.names = TRUE,quote=FALSE)

md<-import_qiime_sample_data("/mnt/genomics/Stephany/liverpath/LCMS/metadata_table/SFR21_0514_metadata_table_updatedTIPS_long_v3.txt")
tree<-read_qza("qtree_analysis/qemistree.qza")$data

OTU = otu_table(spectral_counts, taxa_are_rows = TRUE)
TAX = tax_table(feature_table)
physeq = phyloseq(OTU, TAX, md, tree)

liver_diag <- physeq %>% subset_samples(ATTRIBUTE_blood_procedure %in% c("pre","post","bd","return"))
physeq_class <- liver_diag %>% tax_glom("class")
dv_physeq_class <- divnet(physeq_class,ncores = 8,B=20)
#dv_physeq_class <- divnet(physeq_class,X="ATTRIBUTE_blood_procedure",ncores = 8,B=20)
saveRDS(dv_physeq_class,"diversity_metrics/divnet_HElcms_sampl.rds")
#saveRDS(dv_physeq_class,"diversity_metrics/divnet_HElcms.rds")
dv_physeq_class<-readRDS("diversity_metrics/divnet_HElcms.rds")

write.table(dv_physeq_class$shannon %>%summary,"diversity_metrics/SFR21_0517_divnet_shannonSumm_sampllevel.txt",sep = "\t",row.names = TRUE,quote=FALSE)
write.table(dv_physeq_class$simpson %>%summary,"diversity_metrics/SFR21_0517_divnet_simpsonSumm_sampllevel.txt",sep = "\t",row.names = TRUE,quote=FALSE)
write.table(dv_physeq_class$`bray-curtis`,"diversity_metrics/SFR21_0517_divnet_braycurtisCorr.txt",sep = "\t",row.names = TRUE,quote=FALSE)

p1<-plot(dv_physeq_class$shannon, physeq_class, col = "ATTRIBUTE_blood_procedure")
ggsave("diversity_metrics/SFR21_0517_shannonbysampl_samplevel.png", plot=p1,height=6, width=6)
p2<-plot(dv_physeq_class$simpson, physeq_class, col = "ATTRIBUTE_blood_procedure")
ggsave("diversity_metrics/SFR21_0517_simpsonsampl_samplelevel.png", plot=p2,height=6, width=6)

testDiversity(dv_physeq_class, "shannon")
testDiversity(dv_physeq_class, "simpson")

##################################################
md<-fread("/data/Stephany/liverpath/data_files/LCMS/metadata_table/metadata_table-00000.tsv")%>%
  rename(sample_name=filename)
shannon<-fread("diversity_metrics/SFR21_0517_divnet_shannonSumm_sampllevel.txt")%>%left_join(.,md,by="sample_names")%>%
  mutate(ATTRIBUTE_blood_procedure=factor(ATTRIBUTE_blood_procedure,c("pre","post","bd","return")))

p<-ggplot(shannon, aes(x=ATTRIBUTE_blood_procedure, y=estimate)) + geom_boxplot()+geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  labs(y="shannon estimate",x="TIPS visit")+theme_minimal()
ggsave("diversity_metrics/shannon_boxplotTIPSvisit.png", plot=p,height=3, width=4)
pairwise.wilcox.test(shannon$estimate,shannon$ATTRIBUTE_blood_procedure,p.adjust.method = "BH")

p<-ggplot(shannon, aes(x=blood_draw, y=estimate)) + geom_boxplot()+geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  labs(y="shannon estimate",x="serum source")+theme_minimal()
ggsave("diversity_metrics/shannon_boxplotserum_source.png", plot=p,height=3, width=4)
pairwise.wilcox.test(shannon$estimate,shannon$blood_draw,p.adjust.method = "BH")

p<-ggplot(shannon, aes(x=factor(TIPS_HE), y=estimate)) + geom_boxplot()+geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  labs(y="shannon estimate",x="HE diagnosis")+theme_minimal()
ggsave("diversity_metrics/shannon_boxplotHEdiagnosis.png", plot=p,height=3, width=4)
pairwise.wilcox.test(shannon$estimate,shannon$TIPS_HE,p.adjust.method = "BH")

p<-ggplot(shannon, aes(x=ATTRIBUTE_diagnosis_short, y=estimate)) + geom_boxplot()+geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  labs(y="shannon estimate",x="other conditions")+theme_minimal()+theme(axis.text.x = element_text(angle = 30))
ggsave("diversity_metrics/shannon_boxplotothercond.png", plot=p,height=4, width=7)
pairwise.wilcox.test(shannon$estimate,shannon$ATTRIBUTE_diagnosis_short,p.adjust.method = "BH")

p<-ggplot(shannon, aes(x=ATTRIBUTE_diagnosis_short, y=estimate)) + geom_boxplot()+geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  labs(y="shannon estimate",x="other conditions")+theme_minimal()+theme(axis.text.x = element_text(angle = 30))
ggsave("diversity_metrics/shannon_boxplotothercond.png", plot=p,height=4, width=6)

p<-ggplot(shannon, aes(x=ammonia, y=estimate,color=ATTRIBUTE_blood_procedure)) + geom_point()+theme_minimal()+theme(legend.position = "none")
ggsave("diversity_metrics/shannon_boxplot_ammonia.png", plot=p,height=4, width=6)

##################################################

#plot all the mtbs
mtbid_dict<-fread("/data/Stephany/liverpath/data_files/LCMS/DB_result/32405003deff48dfb9e7cb571ecefc23.tsv")%>%
  rename(rowID=`#Scan#`)

md<-fread("/data/Stephany/liverpath/data_files/LCMS/metadata_table/metadata_table-00000.tsv")%>%rename(sample_name=filename)
# md$ATTRIBUTE_blood_procedure[md$ATTRIBUTE_Sample=="Pre-HET"]<-"pre"
# md$ATTRIBUTE_blood_procedure[grepl("Post-",md$ATTRIBUTE_Sample)]<-"post"
# md$ATTRIBUTE_blood_procedure[grepl("POD-1",md$ATTRIBUTE_Sample)]<-"post"
md$ATTRIBUTE_blood_procedure<-factor(md$ATTRIBUTE_blood_procedure,c("pre","post","bd","return"))


mtb_all<-fread("data_files/LCMS/quantification_table/quantification_table-00000.csv")%>%
  gather(ATTRIBUTE_sampleID,counts, -`row ID`,-`row m/z`,-`row retention time`)%>% rename(rowID=`row ID`)%>%
  mutate(log_counts=log10(counts+1))%>% left_join(.,md, by="ATTRIBUTE_sampleID")%>%
  left_join(.,mtbid_dict,by="rowID")

write.table(mtb_all,"diff_HE_analysis/SFR21_1025_HEmtbdata_allwmd.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#find mtbs with sig differences
mtb_periph<-mtb_all%>%filter(ATTRIBUTE_blood_origin=="peripheral")%>%as.data.table()
mtb_hv<-mtb_all%>%filter(ATTRIBUTE_blood_origin=="hepatic")%>%as.data.table()

run_wilcox<-function (mtb_auc.dt) {
  pr_mtbids.mlms = lapply(mtb_auc.dt[, unique(mtb_id)], function (oid) {
    cat('mtb_ID: ', oid, '\n')
    mtb_auc_sub.dt = mtb_auc.dt[mtb_id == oid]
    if (nrow(mtb_auc_sub.dt) > 3) {
      if(length(unique(mtb_auc_sub.dt$ATTRIBUTE_blood_procedure))>1){
        #tryCatch(wilcox.test(mtb_auc_sub.dt$log_counts ~ mtb_auc_sub.dt$ATTRIBUTE_blood_procedure), error = function(e) NULL)
        #tryCatch(pairwise.wilcox.test(mtb_auc_sub.dt$log_counts, mtb_auc_sub.dt$ATTRIBUTE_blood_procedure,p.adjust.method="bonferroni"), error = function(e) NULL)
        tryCatch(pairwise.wilcox.test(mtb_auc_sub.dt$log_counts, mtb_auc_sub.dt$ATTRIBUTE_blood_procedure,p.adjust.method="fdr"), error = function(e) NULL)
      }
    }
  })
  names(pr_mtbids.mlms) = mtb_auc.dt[, unique(mtb_id)]
  return(pr_mtbids.mlms)
}

getpvals<-function(m_m1){
  m_m1<-plyr::compact(m_m1)
  pr_mtbid_anova_pvals = lapply(m_m1, function(x) x$p.value)
  osuids_sig_pvals=t(matrix(unlist(pr_mtbid_anova_pvals[names(m1)]),nrow=4))
  #osuids_sig_pvals=t(matrix(unlist(pr_mtbid_anova_pvals[names(m1)]),nrow=1))
  colnames(osuids_sig_pvals) <- c("post_vs_bd","pre_vs_bd","post_vs_post","pre_vs_post")
  #colnames(osuids_sig_pvals) <- c("pre_vs_post")
  osuids_sig_pvals<-as.data.table(osuids_sig_pvals)%>%mutate(mtb_id=names(m_m1))%>%select(mtb_id,everything(),-post_vs_post)%>%
    gather(comparison,pval,-mtb_id)
  
  return(osuids_sig_pvals)
}

#HE diag vs. pre/post/bd
m1<-run_wilcox(mtb_PHEdiag)
ppvals_Wilcox<-getpvals(m1)
sigmtbs<-ppvals_Wilcox%>%filter(pval<0.05)
sigmtb_all<-mtb_PHEdiag%>%filter(mtb_id %in% sigmtbs$mtb_id)%>%
  mutate(label_name=paste(mtb_id,direct_parent,sep=" "),
         ATTRIBUTE_blood_procedure=factor(ATTRIBUTE_blood_procedure, levels = c("pre", "post", "bd")),
         mtb_id=as.numeric(mtb_id))
labels<-unique(sigmtb_all$label_name)
names(labels)<-unique(sigmtb_all$mtb_id)
my_comparisons <- list( c("pre", "post"), c("pre", "bd"), c("post", "bd"))
#my_comparisons <- list( c("pre", "post"))
p <- ggplot(sigmtb_all, aes(x=ATTRIBUTE_blood_procedure, y=log_counts)) + geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+ #geom_violin()+geom_boxplot(width=0.1)+#
  facet_wrap(~mtb_id,scale="free",labeller=labeller(mtb_id=labels))+stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_minimal()+labs(x="HE visit", y="log10 mtb abun",title="Hepatic Vein Blood mtbs of Participants that Developed HE")
ggsave("diff_HE_analysis/SFR21_0628_violinplot_HeVisitvslogMtbAbun_HVDevHE_sigwilcoxFDR_pless0.0.05.pdf", plot=p,height=14, width=13)
ggsave("diff_HE_analysis/SFR21_0628_boxplot_HeVisitvslogMtbAbun_HVDevHE_sigwilcoxFDR_pless0.0.05.pdf", plot=p,height=14, width=13)

#HE constant pre/post/bd
mtb_PHEconst<-mtb_periph%>%filter(HE_outcome=="constant_HE"|HE_outcome=="worst_HE")
#mtb_PHEconst<-mtb_hv%>%filter(HE_outcome=="constant_HE"|HE_outcome=="worst_HE")

m2<-run_wilcox(mtb_PHEconst)
ppvals_Wilcox2<-getpvals(m2)
sigmtbs2<-ppvals_Wilcox2%>%filter(pval<0.05)
sigmtb_all<-mtb_PHEconst%>%filter(mtb_id %in% sigmtbs2$mtb_id)%>%
  mutate(label_name=paste(mtb_id,direct_parent,sep=" "),
         ATTRIBUTE_blood_procedure=factor(ATTRIBUTE_blood_procedure, levels = c("pre", "post", "bd")),
         mtb_id=as.numeric(mtb_id))
labels<-unique(sigmtb_all$label_name)
names(labels)<-unique(sigmtb_all$mtb_id)
my_comparisons <- list( c("pre", "post"), c("pre", "bd"), c("post", "bd"))
#my_comparisons <- list( c("pre", "post"))
p <- ggplot(sigmtb_all, aes(x=ATTRIBUTE_blood_procedure, y=log_counts)) +geom_violin()+geom_boxplot(width=0.1)+# geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+#
  facet_wrap(~mtb_id,scale="free",labeller=labeller(mtb_id=labels))+stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_minimal()+labs(x="HE visit", y="log10 mtb abun",title="Hepatic Vein Blood mtbs of Participants with Persistent or Worsening HE")
ggsave("diff_HE_analysis/SFR21_0628_violinplot_HeVisitvslogMtbAbun_HVConstHE_sigwilcoxFDR_pless0.0.05.pdf", plot=p,height=10, width=13)
#ggsave("diff_HE_analysis/SFR21_0628_boxplot_HeVisitvslogMtbAbun_HVConstHE_sigwilcoxFDR_pless0.0.05.pdf", plot=p,height=10, width=13)

#no HE pre/post/bd
mtb_PHEno<-mtb_periph%>%filter(HE_outcome=="no_HE")
#mtb_PHEno<-mtb_hv%>%filter(HE_outcome=="no_HE")
m3<-run_wilcox(mtb_PHEno)
ppvals_Wilcox3<-getpvals(m3)
sigmtbs3<-ppvals_Wilcox3%>%filter(pval<0.05)
sigmtb_all<-mtb_PHEno%>%filter(mtb_id %in% sigmtbs3$mtb_id)%>%
  mutate(label_name=paste(mtb_id,direct_parent,sep=" "),
         ATTRIBUTE_blood_procedure=factor(ATTRIBUTE_blood_procedure, levels = c("pre", "post", "bd")),
         mtb_id=as.numeric(mtb_id))
labels<-unique(sigmtb_all$label_name)
names(labels)<-unique(sigmtb_all$mtb_id)
my_comparisons <- list( c("pre", "post"), c("pre", "bd"), c("post", "bd"))
#my_comparisons <- list( c("pre", "post"))
p <- ggplot(sigmtb_all, aes(x=ATTRIBUTE_blood_procedure, y=log_counts)) +geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+#geom_violin()+geom_boxplot(width=0.1)+
  facet_wrap(~mtb_id,scale="free",labeller=labeller(mtb_id=labels))+stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_minimal()+labs(x="HE visit", y="log10 mtb abun",title="Hepatic Vein Blood mtbs of Participants with no HE Development")
ggsave("diff_HE_analysis/SFR21_0628_violinplot_HeVisitvslogMtbAbun_HVNoHE_sigwilcoxFDR_pless0.0.05.pdf", plot=p,height=8, width=8)
ggsave("diff_HE_analysis/SFR21_0628_boxplot_HeVisitvslogMtbAbun_HVNoHE_sigwilcoxFDR_pless0.0.05.pdf", plot=p,height=8, width=8)

ppvals_Wilcox<-ppvals_Wilcox%>%mutate(He_outcome="Developed_HE")
ppvals_Wilcox2<-ppvals_Wilcox2%>%mutate(He_outcome="Persistent_HE")
ppvals_Wilcox3<-ppvals_Wilcox3%>%mutate(He_outcome="No_HE")
all_ppvals<-rbind(ppvals_Wilcox,ppvals_Wilcox2,ppvals_Wilcox3)

#write.table(all_ppvals,"diff_HE_analysis/SFR21_0625_Perif_HEVisitvslogMtbAbun_wilcoxFDR_allHEoutcomes.txt",sep = "\t",row.names = FALSE,quote=FALSE)
write.table(all_ppvals,"diff_HE_analysis/SFR21_0628_HV_HEVisitvslogMtbAbun_wilcoxFDR_allHEoutcomes.txt",sep = "\t",row.names = FALSE,quote=FALSE)
##################################################

#venn diagram of mtbs with sig differences for HE outcome
sigmtbs<-sigmtbs%>%mutate(He_outcome="Developed_HE")
sigmtbs2<-sigmtbs2%>%mutate(He_outcome="Persistent_HE")
sigmtbs3<-sigmtbs3%>%mutate(He_outcome="no_HE")
sigmtbs_allP<-rbind(sigmtbs,sigmtbs2,sigmtbs3)
#sigmtbs_allhv<-rbind(sigmtbs,sigmtbs2,sigmtbs3)
write.table(sigmtbs_allP,"diff_HE_analysis/SFR21_0628_Perif_HEVisitvslogMtbAbun_sigwilcoxFDRpless0.05_allHEoutcomes.txt",sep = "\t",row.names = FALSE,quote=FALSE)
#write.table(sigmtbs_allhv,"diff_HE_analysis/SFR21_0628_HV_HEVisitvslogMtbAbun_sigwilcoxFDRpless0.05_allHEoutcomes.txt",sep = "\t",row.names = FALSE,quote=FALSE)

prepost1<-(sigmtbs_allhv%>%filter(He_outcome=="Developed_HE"))$mtb_id
prepost2<-(sigmtbs_allhv%>%filter(He_outcome=="Persistent_HE"))$mtb_id
prepost3<-(sigmtbs_allhv%>%filter(He_outcome=="No_HE"))$mtb_id

venn.diagram(
  x = list(prepost1,prepost2,prepost3),
  category.names = c("Developed_HE" , "Persistent_HE" , "No_HE"),
  #filename = 'diff_HE_analysis/SFR21_0608_venn_Perif_HEVisitvslogMtbAbun_sigwilcoxpless0.05_prepost.png',
  filename = 'diff_HE_analysis/SFR21_0628_venn_HV_HEVisitvslogMtbAbun_sigwilcoxpless0.05_prepost.png',
  output=TRUE
)

sigmtbs_allP_prepost<-sigmtbs_allP%>%filter(comparison=="pre_vs_post")
prepost1<-(sigmtbs_allP_prepost%>%filter(He_outcome=="Developed_HE"))$mtb_id
prepost2<-(sigmtbs_allP_prepost%>%filter(He_outcome=="Persistent_HE"))$mtb_id
prepost3<-(sigmtbs_allP_prepost%>%filter(He_outcome=="no_HE"))$mtb_id

venn.diagram(
  x = list(prepost1,prepost2,prepost3),
  category.names = c("Developed_HE" , "Persistent_HE" , "no_HE"),
  filename = 'diff_HE_analysis/SFR21_0628_venn_Perif_HEVisitvslogMtbAbun_sigwilcoxFDRpless0.05_prepost.png',
  output=TRUE
)

sigmtbs_allP_prebd<-sigmtbs_allP%>%filter(comparison=="pre_vs_bd")
prebd1<-(sigmtbs_allP_prebd%>%filter(He_outcome=="Developed_HE"))$mtb_id
prebd2<-(sigmtbs_allP_prebd%>%filter(He_outcome=="Persistent_HE"))$mtb_id
prebd3<-(sigmtbs_allP_prebd%>%filter(He_outcome=="no_HE"))$mtb_id

venn.diagram(
  x = list(prebd1,prebd2,prebd3),
  category.names = c("Developed_HE" , "Persistent_HE" , "no_HE"),
  filename = 'diff_HE_analysis/SFR21_0628_venn_Perif_HEVisitvslogMtbAbun_sigwilcoxFDRpless0.05_prebd.png',
  output=TRUE
)

sigmtbs_allP_postbd<-sigmtbs_allP%>%filter(comparison=="post_vs_bd")
postbd1<-(sigmtbs_allP_postbd%>%filter(He_outcome=="Developed_HE"))$mtb_id
postbd2<-(sigmtbs_allP_postbd%>%filter(He_outcome=="Persistent_HE"))$mtb_id
postbd3<-(sigmtbs_allP_postbd%>%filter(He_outcome=="no_HE"))$mtb_id

venn.diagram(
  x = list(postbd1,postbd2,postbd3),
  category.names = c("Developed_HE" , "Persistent_HE" , "no_HE"),
  filename = 'diff_HE_analysis/SFR21_0628_venn_Perif_HEVisitvslogMtbAbun_sigwilcoxFDRpless0.05_postbd.png',
  output=TRUE
)

##################################################
sigmtbs_allP<-fread("diff_HE_analysis/SFR21_0628_Perif_HEVisitvslogMtbAbun_sigwilcoxFDRpless0.05_allHEoutcomes.txt")
sigmtbs_allhv<-fread("diff_HE_analysis/SFR21_0628_HV_HEVisitvslogMtbAbun_sigwilcoxFDRpless0.05_allHEoutcomes.txt")

prepost<-sigmtbs_allhv%>%filter(comparison=="pre_vs_post")
y<-list((prepost%>%filter(He_outcome=="Developed_HE"))$mtb_id,
(prepost%>%filter(He_outcome=="Persistent_HE"))$mtb_id,(prepost%>%filter(He_outcome=="No_HE"))$mtb_id)

# prebd<-sigmtbs_allP%>%filter(comparison=="pre_vs_bd")
# y<-list((prebd%>%filter(He_outcome=="Developed_HE"))$mtb_id,
#         (prebd%>%filter(He_outcome=="Persistent_HE"))$mtb_id,(prebd%>%filter(He_outcome=="no_HE"))$mtb_id)

# library(RVenn)
 x<-Venn(y)
 overlap_pairs(x)
 overlap(x)
 discern(x,"Set_1","Set_2")

mtb_all<-fread("diff_HE_analysis/SFR21_0608_HEmtbdata_allwmd.txt")
mtb_periph<-mtb_all%>%filter(blood_draw=="P")%>%
  filter(HE_outcome=="developed_HE_1"|HE_outcome=="developed_HE_2")%>%as.data.table()
mtb_hv<-mtb_all%>%filter(blood_draw=="HV")%>%filter(HE_outcome=="developed_HE_1"|HE_outcome=="developed_HE_2")%>%as.data.table()

mtb_ints<-c(146,229,261,52,13,36,39,110,154,184,185,225,228,234,257,72)
mtb_ints<-c(92,115,117,146,153,171,201,229,261,256,8,52,250,262)
#mtb_ints<-c(115,117,146,153,171,256,8,52,262)
#mtb_ints<-c(13,72,146,228,52,36,39,110,154,184,185,225,229,234,257,261)
mtb_ints_df<-mtb_hv%>%filter(mtb_id %in% mtb_ints)%>%
  mutate(ATTRIBUTE_blood_procedure=factor(ATTRIBUTE_blood_procedure, levels = c("pre", "post", "bd")),
         mtb_id=factor(mtb_id, levels = c(146,229,261,52,13,36,39,110,154,184,185,225,228,234,257,72)),
         #mtb_id=factor(mtb_id, levels = c(92,115,117,146,153,171,201,229,261,256,8,52,250,262)),
         label_name=paste(mtb_id,direct_parent,sep=" "))

#mtb_ints_df$GNPS_match[mtb_ints_df$mtb_id==115]<-"1,2-dioleoyl-sn-glycero-3-phosphatidylcholine"
# mtb_ints_df$GNPS_match[mtb_ints_df$mtb_id==228]<-"Sulfadoxine"
# mtb_ints_df$GNPS_match[mtb_ints_df$mtb_id==52]<-"(2-Chloro-5-nitrophenyl)-[4-(2-nitrophenyl)sulfonylpiperazin-1-yl]methanone"
# mtb_ints_df$GNPS_match[mtb_ints_df$mtb_id==36]<-"Diethyl 2-(2,2-diethoxyethenyl)butanedioate"
# mtb_ints_df$GNPS_match[mtb_ints_df$mtb_id==185]<-"3,7-Dimethyl-5,9-dioxo-11-hydroxy-4,8-dioxadodecanoic acid ethyl ester"
# mtb_ints_df$GNPS_match[mtb_ints_df$mtb_id==257]<-"3-Hexadecyl-3-propylnonadecan-1-amine"
labels<-unique(mtb_ints_df$label_name)
names(labels)<-unique(mtb_ints_df$mtb_id)
#my_comparisons <- list( c("pre", "post"), c("pre", "bd"), c("post", "bd"))
my_comparisons <- list( c("pre", "post"))
p <- ggplot(mtb_ints_df, aes(x=ATTRIBUTE_blood_procedure, y=log_counts)) +geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+ # geom_violin()+geom_boxplot(width=0.1)+
  facet_wrap(~mtb_id,scale="free",ncol=5,labeller=labeller(mtb_id=labels))+stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_minimal()+labs(x="HE visit", y="log10 mtb abun",title="Interesting Peripheral mtbs of Participants that Developed HE")
ggsave("diff_HE_analysis/SFR21_0628_boxplot_HeVisitvslogMtbAbun_PerifDevHE_sigwilcoxFDR_pless0.0.05_uniqueIntmtbs.pdf", plot=p,height=11, width=12.5)

p <- ggplot(mtb_ints_df, aes(x=ATTRIBUTE_blood_procedure, y=log_counts)) +geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+ # geom_violin()+geom_boxplot(width=0.1)+
  facet_wrap(~mtb_id,scale="free",ncol=5,labeller=labeller(mtb_id=labels))+stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_minimal()+labs(x="HE visit", y="log10 mtb abun",title="Interesting Hepatic Vein mtbs of Participants that Developed HE")
ggsave("diff_HE_analysis/SFR21_0628_boxplot_HeVisitvslogMtbAbun_HVDevHE_sigwilcoxFDR_pless0.0.05_uniqueIntmtbs.pdf", plot=p,height=12, width=13)
##################################################
#combining peripheral with HV
mtb_all<-fread("diff_HE_analysis/SFR21_0608_HEmtbdata_allwmd.txt")
run_wilcox<-function (mtb_auc.dt) {
  pr_mtbids.mlms = lapply(mtb_auc.dt[, unique(mtb_id)], function (oid) {
    cat('mtb_ID: ', oid, '\n')
    mtb_auc_sub.dt = mtb_auc.dt[mtb_id == oid]
    if (nrow(mtb_auc_sub.dt) > 3) {
      if(length(unique(mtb_auc_sub.dt$ATTRIBUTE_blood_procedure))>1){
        tryCatch(pairwise.wilcox.test(mtb_auc_sub.dt$log_counts, mtb_auc_sub.dt$ATTRIBUTE_blood_procedure,p.adjust.method="bonferroni"), error = function(e) NULL)
        #tryCatch(pairwise.wilcox.test(mtb_auc_sub.dt$log_counts, mtb_auc_sub.dt$ATTRIBUTE_blood_procedure,p.adjust.method="fdr"), error = function(e) NULL)
      }
    }
  })
  names(pr_mtbids.mlms) = mtb_auc.dt[, unique(mtb_id)]
  return(pr_mtbids.mlms)
}

getpvals<-function(m_m1){
  m_m1<-plyr::compact(m_m1)
  pr_mtbid_anova_pvals = lapply(m_m1, function(x) x$p.value)
  osuids_sig_pvals=t(matrix(unlist(pr_mtbid_anova_pvals[names(m1)]),nrow=4))
  colnames(osuids_sig_pvals) <- c("post_vs_bd","pre_vs_bd","post_vs_post","pre_vs_post")
  osuids_sig_pvals<-as.data.table(osuids_sig_pvals)%>%mutate(mtb_id=names(m_m1))%>%select(mtb_id,everything(),-post_vs_post)%>%
    gather(comparison,pval,-mtb_id)
  
  return(osuids_sig_pvals)
}

#dev HE
mtb_PHEdiag<-mtb_all%>%filter(HE_outcome=="developed_HE"| HE_outcome=="developed_HE_1"|HE_outcome=="developed_HE_2")
m1<-run_wicox(mtb_PHEdiag)
ppvals_Wilcox<-getpvals(m1)
sigmtbs<-ppvals_Wilcox%>%filter(pval<0.05)
sigmtb_all<-mtb_PHEdiag%>%filter(mtb_id %in% sigmtbs$mtb_id)%>%
  mutate(label_name=paste(mtb_id,direct_parent,sep=" "),
         ATTRIBUTE_blood_procedure=factor(ATTRIBUTE_blood_procedure, levels = c("pre", "post", "bd")),
         mtb_id=as.numeric(mtb_id))
labels<-unique(sigmtb_all$label_name)
names(labels)<-unique(sigmtb_all$mtb_id)
my_comparisons <- list( c("pre", "post"), c("pre", "bd"), c("post", "bd"))
p <- ggplot(sigmtb_all, aes(x=ATTRIBUTE_blood_procedure, y=log_counts)) + geom_violin()+geom_boxplot(width=0.1)+#geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+ #
  facet_wrap(~mtb_id,scale="free",labeller=labeller(mtb_id=labels))+stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_minimal()+labs(x="HE visit", y="log10 mtb abun",title="mtbs of Participants that Developed HE")
ggsave("diff_HE_analysis/SFR21_0629_violinplot_HeVisitvslogMtbAbun_HVDevHE_sigwilcoxBonf_pless0.0.05.pdf", plot=p,height=20, width=16)
ggsave("diff_HE_analysis/SFR21_0629_boxplot_HeVisitvslogMtbAbun_HVDevHE_sigwilcoxBonf_pless0.0.05.pdf", plot=p,height=20, width=16)

#Persist HE
mtb_PHEconst<-mtb_all%>%filter(HE_outcome=="constant_HE"|HE_outcome=="worst_HE")
m2<-run_wilcox(mtb_PHEconst)
ppvals_Wilcox2<-getpvals(m2)
sigmtbs2<-ppvals_Wilcox2%>%filter(pval<0.05)
sigmtb_all<-mtb_PHEconst%>%filter(mtb_id %in% sigmtbs2$mtb_id)%>%
  mutate(label_name=paste(mtb_id,direct_parent,sep=" "),
         ATTRIBUTE_blood_procedure=factor(ATTRIBUTE_blood_procedure, levels = c("pre", "post", "bd")),
         mtb_id=as.numeric(mtb_id))
labels<-unique(sigmtb_all$label_name)
names(labels)<-unique(sigmtb_all$mtb_id)
my_comparisons <- list( c("pre", "post"), c("pre", "bd"), c("post", "bd"))
p <- ggplot(sigmtb_all, aes(x=ATTRIBUTE_blood_procedure, y=log_counts)) +geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+#geom_violin()+geom_boxplot(width=0.1)+# 
  facet_wrap(~mtb_id,scale="free",labeller=labeller(mtb_id=labels))+stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_minimal()+labs(x="HE visit", y="log10 mtb abun",title="mtbs of Participants with Persistent or Worsening HE")
ggsave("diff_HE_analysis/SFR21_0629_violinplot_HeVisitvslogMtbAbun_ConstHE_sigwilcoxBonf_pless0.0.05.pdf", plot=p,height=20, width=16)
ggsave("diff_HE_analysis/SFR21_0629_boxplot_HeVisitvslogMtbAbun_ConstHE_sigwilcoxBonf_pless0.0.05.pdf", plot=p,height=20, width=16)

#no HE
mtb_PHEno<-mtb_all%>%filter(HE_outcome=="no_HE")
m3<-run_wilcox(mtb_PHEno)
ppvals_Wilcox3<-getpvals(m3)
sigmtbs3<-ppvals_Wilcox3%>%filter(pval<0.05)
sigmtb_all<-mtb_PHEno%>%filter(mtb_id %in% sigmtbs3$mtb_id)%>%
  mutate(label_name=paste(mtb_id,direct_parent,sep=" "),
         ATTRIBUTE_blood_procedure=factor(ATTRIBUTE_blood_procedure, levels = c("pre", "post", "bd")),
         mtb_id=as.numeric(mtb_id))
labels<-unique(sigmtb_all$label_name)
names(labels)<-unique(sigmtb_all$mtb_id)
my_comparisons <- list( c("pre", "post"), c("pre", "bd"), c("post", "bd"))
p <- ggplot(sigmtb_all, aes(x=ATTRIBUTE_blood_procedure, y=log_counts)) +geom_violin()+geom_boxplot(width=0.1)+#geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+#
  facet_wrap(~mtb_id,scale="free",labeller=labeller(mtb_id=labels))+stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme_minimal()+labs(x="HE visit", y="log10 mtb abun",title="mtbs of Participants with no HE Development")
ggsave("diff_HE_analysis/SFR21_0629_violinplot_HeVisitvslogMtbAbun_NoHE_sigwilcoxBonf_pless0.0.05.pdf", plot=p,height=16, width=14)
ggsave("diff_HE_analysis/SFR21_0629_boxplot_HeVisitvslogMtbAbun_NoHE_sigwilcoxBonf_pless0.0.05.pdf", plot=p,height=16, width=14)

##################################################
mtb_all<-fread("diff_HE_analysis/SFR21_0608_HEmtbdata_allwmd.txt")
BAs<-fread("diff_HE_analysis/BA_network_IDs.txt")

# pairedsamps<-mtb_pp%>%group_by(mtb_id,ATTRIBUTE_unique_pt_id,blood_draw,ATTRIBUTE_blood_procedure)%>%summarise(n=n())%>%
#   filter(blood_draw=="P"|blood_draw=="HV")%>%spread(ATTRIBUTE_blood_procedure,n)%>%
#   group_by(mtb_id,blood_draw)%>%summarise(pre_n=sum(pre,na.rm=TRUE),post_n=sum(post,na.rm=TRUE))%>%
#   mutate(remove=pre_n-post_n)
# mtb_pp_prd<-mtb_pp%>%left_join(.,pairedsamps, by=c("mtb_id","ATTRIBUTE_unique_pt_id","blood_draw"))%>%filter(paired=="yes")
  

##################################################
#looking at distance metric diffs based on HE grading (within indv)

#md<-fread("data_files/LCMS/metadata_table/metadata_table-peripheral-prepost.txt") %>%
md<-fread("data_files/LCMS/metadata_table/metadata_table-hepatic.txt") %>%
#md<-fread("data_files/LCMS/metadata_table/metadata_table-hepatic-rmoutlier.txt") %>%
#md<-fread("data_files/LCMS/metadata_table/metadata_table-peripheral.txt") %>%
  #filter(ATTRIBUTE_blood_procedure!="post")%>%
  dplyr::select(sampleid,ATTRIBUTE_pt_id, ATTRIBUTE_blood_procedure, ATTRIBUTE_PreTIPS_HE_mod,ATTRIBUTE_Worst_PostTIPS_HE_mod)%>%
  spread(ATTRIBUTE_blood_procedure,sampleid)

#perif_dist<-fread("diversity_metrics_LCMS/beta_deicode_lcms_peripheral-prepost/distance-matrix.tsv")%>%
perif_dist<-fread("diversity_metrics_LCMS/beta_deicode_lcms_hepatic/distance-matrix.tsv")%>%
#perif_dist<-fread("diversity_metrics_LCMS/beta_deicode_lcms_hepatic-rmoutlier/distance-matrix.tsv")%>%
#perif_dist<-fread("diversity_metrics_LCMS/beta_deicode_lcms_peripheral/distance-matrix.tsv")%>%
  gather(V2,distance,-V1) %>%dplyr::rename(pre=V1, post=V2) %>%
  right_join(.,md,by=c("pre","post")) %>%
  drop_na() %>%
  mutate(ATTRIBUTE_PreTIPS_HE_mod=as.factor(ATTRIBUTE_PreTIPS_HE_mod),
         ATTRIBUTE_Worst_PostTIPS_HE_mod=as.factor(ATTRIBUTE_Worst_PostTIPS_HE_mod))

p<-ggplot(perif_dist, aes(x=ATTRIBUTE_PreTIPS_HE_mod, y=distance,fill=ATTRIBUTE_PreTIPS_HE_mod)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#39B54A", "#283891", "#EF3E36"))+scale_fill_manual(values=c("#39B54A", "#283891", "#EF3E36"))+
  labs(x="HE grade",y="robust aitchison distance", title="Pre TIPS grading")+
  theme(legend.position = "none")

#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_peripheral/distancevs.preTIPS_prebd.pdf", plot=p,height=3, width=2.8)
#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_hepatic/distancevs.preTIPS.pdf", plot=p,height=3, width=2.8)
ggsave("diversity_metrics_LCMS/beta_deicode_lcms_peripheral-prepost/distancevs.preTIPS.pdf", plot=p,height=3, width=2.8)

pairwise.wilcox.test(perif_dist$distance, perif_dist$ATTRIBUTE_PreTIPS_HE_mod,
                     p.adjust.method="fdr")

#perif
# 0   
# 1 0.55

#hep
# 0
# 1 1

#prevsbd
# 0  
# 1 0.6


p<-ggplot(perif_dist, aes(x=fct_reorder(ATTRIBUTE_pt_id,distance), y=distance, fill=ATTRIBUTE_PreTIPS_HE_mod)) +
  geom_bar(stat="identity")+theme_classic() + scale_color_manual(values=c("#39B54A", "#283891", "#EF3E36"))+scale_fill_manual(values=c("#39B54A", "#283891", "#EF3E36"))+
  labs(x="subject ID",y="robust aitchison distance", title="Pre TIPS grading", fill="HE grading")+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_peripheral/distanceprebd_bysubject.pdf", plot=p,height=3.5, width=5)
#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_hepatic/distancepre_bysubject.pdf", plot=p,height=3.5, width=5)
ggsave("diversity_metrics_LCMS/beta_deicode_lcms_peripheral-prepost/distancepre_bysubject.pdf", plot=p,height=3.5, width=5)

p<-ggplot(perif_dist, aes(x=ATTRIBUTE_Worst_PostTIPS_HE_mod, y=distance,fill=ATTRIBUTE_Worst_PostTIPS_HE_mod)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#39B54A", "#283891", "#EF3E36"))+scale_fill_manual(values=c("#39B54A", "#283891", "#EF3E36"))+
  labs(x="HE grade",y="robust aitchison distance", title="Post TIPS grading")+
  theme(legend.position = "none")

#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_peripheral/distancevsprebd.postTIPS.pdf", plot=p,height=3, width=3.5)
#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_hepatic/distancevs.postTIPS.pdf", plot=p,height=3, width=3.5)
ggsave("diversity_metrics_LCMS/beta_deicode_lcms_hepatic-rmoutlier/distancevs.postTIPS.pdf", plot=p,height=3, width=3.5)
#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_peripheral-prepost/distancevs.postTIPS.pdf", plot=p,height=3, width=3.5)

pairwise.wilcox.test(perif_dist$distance, perif_dist$ATTRIBUTE_Worst_PostTIPS_HE_mod,
                     p.adjust.method="fdr")

#perif
# 0    1   
# 1  0.87 -   
#   2+ 0.87 0.87
#hep
# 0     1    
# 1  0.027 -    
#   2+ 0.036 0.272

#hep removing outlier
# 0     1    
# 1  0.036 -    
#   2+ 0.036 0.491

#prebd
# 0    1   
# 1  0.42 -   
#   2+ 0.42 0.87urro

x<-perif_dist%>%filter(ATTRIBUTE_Worst_PostTIPS_HE_mod=="0")
mean(x$distance) #45% diff with 20 patient
sd(x$distance)
y<-perif_dist%>%filter(ATTRIBUTE_Worst_PostTIPS_HE_mod=="1")
mean(y$distance)
sd(y$distance)
z<-perif_dist%>%filter(ATTRIBUTE_Worst_PostTIPS_HE_mod=="2+")
mean(z$distance) #82% difference you can detect with 6 patients; 5 fold difference
sd(z$distance)

p<-ggplot(perif_dist, aes(x=fct_reorder(ATTRIBUTE_pt_id,distance), y=distance, fill=ATTRIBUTE_Worst_PostTIPS_HE_mod)) +
  geom_bar(stat="identity")+theme_classic() + scale_color_manual(values=c("#39B54A", "#283891", "#EF3E36"))+scale_fill_manual(values=c("#39B54A", "#283891", "#EF3E36"))+
  labs(x="subject ID",y="robust aitchison distance", title="Post TIPS grading", fill="HE grading")+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_peripheral/distancepost_prebdbysubject.pdf", plot=p,height=3.5, width=5)
#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_peripheral-prepost/distancepost_bysubject.pdf", plot=p,height=3.5, width=5)
#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_hepatic/distancepost_bysubject.pdf", plot=p,height=3.5, width=5)
ggsave("diversity_metrics_LCMS/beta_deicode_lcms_hepatic-rmoutlier/distancepost_bysubject.pdf", plot=p,height=3.5, width=5)

#look at the delta PC to see if there is separation there
#md<-fread("data_files/LCMS/metadata_table/metadata_table-peripheral-prepost.txt") %>%
#md<-fread("data_files/LCMS/metadata_table/metadata_table-hepatic.txt") %>%
md<-fread("data_files/LCMS/metadata_table/metadata_table-peripheral.txt") %>%
  filter(ATTRIBUTE_blood_procedure!="post")%>%
  select(sampleid,ATTRIBUTE_pt_id, ATTRIBUTE_blood_procedure, ATTRIBUTE_PreTIPS_HE_mod,ATTRIBUTE_Worst_PostTIPS_HE_mod)

#md2<-fread("data_files/LCMS/metadata_table/metadata_table-peripheral-prepost.txt") %>%
#md2<-fread("data_files/LCMS/metadata_table/metadata_table-hepatic.txt") %>%
md2<-fread("data_files/LCMS/metadata_table/metadata_table-peripheral.txt") %>%
  filter(ATTRIBUTE_blood_procedure!="post")%>%
  select(ATTRIBUTE_pt_id, ATTRIBUTE_PreTIPS_HE_mod,ATTRIBUTE_Worst_PostTIPS_HE_mod) %>%unique()

#perif_ord <- read_qza("diversity_metrics_LCMS/ordination-peripheral-prepost.qza")
#perif_ord <- read_qza("diversity_metrics_LCMS/ordination-hepatic.qza")
perif_ord <- read_qza("diversity_metrics_LCMS/ordination-peripheral.qza")

rpca<-perif_ord$data$Vectors %>%
  select(SampleID, PC1) %>%
  filter(SampleID %in% md$sampleid) %>%
  dplyr::rename(sampleid="SampleID")%>%
  left_join(.,md,by="sampleid")%>%
  pivot_wider(id_cols = ATTRIBUTE_pt_id, 
              names_from = ATTRIBUTE_blood_procedure, 
              values_from = c("sampleid", "PC1"))%>%
  drop_na() %>%
  mutate(delta_PC=PC1_bd-PC1_pre)%>%
  left_join(.,md2,by="ATTRIBUTE_pt_id")

p<-ggplot(rpca, aes(x=ATTRIBUTE_Worst_PostTIPS_HE_mod, y=delta_PC,fill=ATTRIBUTE_Worst_PostTIPS_HE_mod)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#DC4D40","#754B27","#60B358"))+scale_fill_manual(values=c("#DC4D40","#754B27","#60B358"))+
  labs(x="HE grade",y=expression(paste(Delta,"PC1")), title="Post TIPS grading")+
  theme(legend.position = "none")
ggsave("diversity_metrics_LCMS/beta_deicode_lcms_peripheral/delPC1vs.prebdpostTIPS.pdf", plot=p,height=3, width=3.5)
#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_peripheral-prepost/delPC1vs.postTIPS.pdf", plot=p,height=3, width=3.5)
#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_hepatic/delPC1vs.postTIPS.pdf", plot=p,height=3, width=3.5)

pairwise.wilcox.test(rpca$delta_PC, rpca$ATTRIBUTE_Worst_PostTIPS_HE_mod,
                     p.adjust.method="fdr")
#perif
# 0    1   
# 1  0.33 -   
#   2+ 0.54 0.54

#hep
# 0    1   
# 1  0.86 -   
#   2+ 0.86 0.86

#prebd
# 0    1   
# 1  0.93 -   
#   2+ 0.93 0.93

p<-ggplot(rpca, aes(x=fct_reorder(ATTRIBUTE_pt_id,delta_PC), y=delta_PC, fill=ATTRIBUTE_Worst_PostTIPS_HE_mod)) +
  geom_bar(stat="identity")+theme_classic() + scale_color_manual(values=c("#DC4D40","#754B27","#60B358"))+scale_fill_manual(values=c("#DC4D40","#754B27","#60B358"))+
  labs(x="subject ID",y=expression(paste(Delta,"PC1")), title="Post TIPS grading", fill="HE grading")+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("diversity_metrics_LCMS/beta_deicode_lcms_peripheral/delPC1post_prebdbysubject.pdf", plot=p,height=3.5, width=5)
#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_peripheral-prepost/delPC1post_bysubject.pdf", plot=p,height=3.5, width=5)
#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_hepatic/delPC1post_bysubject.pdf", plot=p,height=3.5, width=5)

p<-ggplot(rpca, aes(x=as.factor(ATTRIBUTE_PreTIPS_HE_mod), y=delta_PC,fill=as.factor(ATTRIBUTE_PreTIPS_HE_mod))) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#DC4D40","#754B27","#60B358"))+scale_fill_manual(values=c("#DC4D40","#754B27","#60B358"))+
  labs(x="HE grade",y=expression(paste(Delta,"PC1")), title="Pre TIPS grading")+
  theme(legend.position = "none")
ggsave("diversity_metrics_LCMS/beta_deicode_lcms_peripheral/delPC1vs.prebdpreTIPS.pdf", plot=p,height=3, width=3.5)
#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_peripheral-prepost/delPC1vs.preTIPS.pdf", plot=p,height=3, width=3.5)
#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_hepatic/delPC1vs.preTIPS.pdf", plot=p,height=3, width=3.5)

pairwise.wilcox.test(rpca$delta_PC, rpca$ATTRIBUTE_PreTIPS_HE_mod,
                     p.adjust.method="fdr")
#perif
# 0    
# 1 0.095

#hep
# 0   
# 1 0.55

#prebd
# 0   
# 1 0.13

p<-ggplot(rpca, aes(x=fct_reorder(ATTRIBUTE_pt_id,delta_PC), y=delta_PC, fill=as.factor(ATTRIBUTE_PreTIPS_HE_mod))) +
  geom_bar(stat="identity")+theme_classic() + scale_color_manual(values=c("#DC4D40","#754B27","#60B358"))+scale_fill_manual(values=c("#DC4D40","#754B27","#60B358"))+
  labs(x="subject ID",y=expression(paste(Delta,"PC1")), title="Pre TIPS grading", fill="HE grading")+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("diversity_metrics_LCMS/beta_deicode_lcms_peripheral/delPC1pre_prebdbysubject.pdf", plot=p,height=3.5, width=5)
#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_peripheral-prepost/delPC1pre_bysubject.pdf", plot=p,height=3.5, width=5)
#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_hepatic/delPC1pre_bysubject.pdf", plot=p,height=3.5, width=5)

#looking at distance metric diffs based on HE grading (bw indv 0_0 vs. 1_1 and 2_2) Post first

#md<-fread("data_files/LCMS/metadata_table/metadata_table-peripheral-prepost.txt") %>%
md<-fread("data_files/LCMS/metadata_table/metadata_table-hepatic.txt") %>%
  select(sampleid,ATTRIBUTE_pt_id, ATTRIBUTE_blood_procedure, ATTRIBUTE_PreTIPS_HE_mod,ATTRIBUTE_Worst_PostTIPS_HE_mod)

#md1<-fread("data_files/LCMS/metadata_table/metadata_table-peripheral-prepost.txt") %>%
md1<-fread("data_files/LCMS/metadata_table/metadata_table-hepatic.txt") %>%
  select(sampleid,ATTRIBUTE_pt_id, ATTRIBUTE_blood_procedure, ATTRIBUTE_PreTIPS_HE_mod,ATTRIBUTE_Worst_PostTIPS_HE_mod)%>%
  filter(ATTRIBUTE_blood_procedure=="post")%>%
  rename(V1=sampleid, V1_ATTRIBUTE_pt_id=ATTRIBUTE_pt_id,V1_ATTRIBUTE_PreTIPS_HE_mod=ATTRIBUTE_PreTIPS_HE_mod, V1_ATTRIBUTE_Worst_PostTIPS_HE_mod=ATTRIBUTE_Worst_PostTIPS_HE_mod)%>%
  select(-ATTRIBUTE_blood_procedure)

#md2<-fread("data_files/LCMS/metadata_table/metadata_table-peripheral-prepost.txt") %>%
md2<-fread("data_files/LCMS/metadata_table/metadata_table-hepatic.txt") %>%
  select(sampleid,ATTRIBUTE_pt_id, ATTRIBUTE_blood_procedure, ATTRIBUTE_PreTIPS_HE_mod,ATTRIBUTE_Worst_PostTIPS_HE_mod)%>%
  filter(ATTRIBUTE_blood_procedure=="post")%>%
  rename(V2=sampleid, V2_ATTRIBUTE_pt_id=ATTRIBUTE_pt_id,V2_ATTRIBUTE_PreTIPS_HE_mod=ATTRIBUTE_PreTIPS_HE_mod, V2_ATTRIBUTE_Worst_PostTIPS_HE_mod=ATTRIBUTE_Worst_PostTIPS_HE_mod)%>%
  select(-ATTRIBUTE_blood_procedure)

postIds<-(md%>%filter(ATTRIBUTE_blood_procedure=="post"))$sampleid

#perif_dist_po<-fread("diversity_metrics_LCMS/beta_deicode_lcms_peripheral-prepost/distance-matrix.tsv")%>%
perif_dist_po<-fread("diversity_metrics_LCMS/beta_deicode_lcms_hepatic/distance-matrix.tsv")%>%
  gather(V2,distance,-V1) %>% 
  filter(V1 %in% postIds) %>%
  filter(V2 %in% postIds) %>%
  left_join(.,md1, by="V1") %>%
  left_join(.,md2, by="V2") %>%
  filter(!(V1_ATTRIBUTE_pt_id==V2_ATTRIBUTE_pt_id)) %>%
  mutate(HE_Wpost_comb=paste(V1_ATTRIBUTE_Worst_PostTIPS_HE_mod,V2_ATTRIBUTE_Worst_PostTIPS_HE_mod, sep="_")) %>%
  group_by(distance) %>% 
  filter(row_number()==1) %>%
  filter(HE_Wpost_comb %in% c("0_0","1_1","2+_2+"))

p<-ggplot(perif_dist_po, aes(x=HE_Wpost_comb, y=distance,fill=HE_Wpost_comb)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#DC4D40","#754B27","#60B358"))+scale_fill_manual(values=c("#DC4D40","#754B27","#60B358"))+
  labs(x="HE grade",y="between subject robust aitchison", title="Post TIPS")+
  theme(legend.position = "none")
#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_peripheral-prepost/postdistancevs.postTIPS_bwsubj.pdf", plot=p,height=3, width=3.5)
ggsave("diversity_metrics_LCMS/beta_deicode_lcms_hepatic/postdistancevs.postTIPS_bwsubj.pdf", plot=p,height=3, width=3.5)

pairwise.wilcox.test(perif_dist_po$distance, perif_dist_po$HE_Wpost_comb,
                     p.adjust.method="fdr")

#perif
# 0_0   1_1    
# 1_1   0.314 -      
#   2+_2+ 0.039 9.4e-06

#hep
# 0_0     1_1 
# 1_1   5.5e-05 -   
#   2+_2+ 5.5e-05 0.1

#looking at distance metric diffs based on HE grading (bw indv 0_0 vs. 1_1 and 2_2) now pre

#md<-fread("data_files/LCMS/metadata_table/metadata_table-peripheral-prepost.txt") %>%
md<-fread("data_files/LCMS/metadata_table/metadata_table-hepatic.txt") %>%
  select(sampleid,ATTRIBUTE_pt_id, ATTRIBUTE_blood_procedure, ATTRIBUTE_PreTIPS_HE_mod,ATTRIBUTE_Worst_PostTIPS_HE_mod)

#md1<-fread("data_files/LCMS/metadata_table/metadata_table-peripheral-prepost.txt") %>%
md1<-fread("data_files/LCMS/metadata_table/metadata_table-hepatic.txt") %>%
  select(sampleid,ATTRIBUTE_pt_id, ATTRIBUTE_blood_procedure, ATTRIBUTE_PreTIPS_HE_mod,ATTRIBUTE_Worst_PostTIPS_HE_mod)%>%
  filter(ATTRIBUTE_blood_procedure=="pre")%>%
  rename(V1=sampleid, V1_ATTRIBUTE_pt_id=ATTRIBUTE_pt_id,V1_ATTRIBUTE_PreTIPS_HE_mod=ATTRIBUTE_PreTIPS_HE_mod, V1_ATTRIBUTE_Worst_PostTIPS_HE_mod=ATTRIBUTE_Worst_PostTIPS_HE_mod)%>%
  select(-ATTRIBUTE_blood_procedure)

#md2<-fread("data_files/LCMS/metadata_table/metadata_table-peripheral-prepost.txt") %>%
md2<-fread("data_files/LCMS/metadata_table/metadata_table-hepatic.txt") %>%
  select(sampleid,ATTRIBUTE_pt_id, ATTRIBUTE_blood_procedure, ATTRIBUTE_PreTIPS_HE_mod,ATTRIBUTE_Worst_PostTIPS_HE_mod)%>%
  filter(ATTRIBUTE_blood_procedure=="pre")%>%
  rename(V2=sampleid, V2_ATTRIBUTE_pt_id=ATTRIBUTE_pt_id,V2_ATTRIBUTE_PreTIPS_HE_mod=ATTRIBUTE_PreTIPS_HE_mod, V2_ATTRIBUTE_Worst_PostTIPS_HE_mod=ATTRIBUTE_Worst_PostTIPS_HE_mod)%>%
  select(-ATTRIBUTE_blood_procedure)

preIds<-(md%>%filter(ATTRIBUTE_blood_procedure=="pre"))$sampleid

#perif_dist_pr<-fread("diversity_metrics_LCMS/beta_deicode_lcms_peripheral-prepost/distance-matrix.tsv")%>%
perif_dist_pr<-fread("diversity_metrics_LCMS/beta_deicode_lcms_hepatic/distance-matrix.tsv")%>%
  gather(V2,distance,-V1) %>% 
  filter(V1 %in% preIds) %>%
  filter(V2 %in% preIds) %>%
  left_join(.,md1, by="V1") %>%
  left_join(.,md2, by="V2") %>%
  filter(!(V1_ATTRIBUTE_pt_id==V2_ATTRIBUTE_pt_id)) %>%
  mutate(HE_pre_comb=paste(V1_ATTRIBUTE_PreTIPS_HE_mod,V1_ATTRIBUTE_PreTIPS_HE_mod, sep="_")) %>%
  group_by(distance) %>% 
  filter(row_number()==1) %>%
  filter(HE_pre_comb %in% c("0_0","1_1"))

p<-ggplot(perif_dist_pr, aes(x=HE_pre_comb, y=distance,fill=HE_pre_comb)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#DC4D40","#754B27","#60B358"))+scale_fill_manual(values=c("#DC4D40","#754B27","#60B358"))+
  labs(x="HE grade",y="between subject robust aitchison", title="Pre TIPS grading")+
  theme(legend.position = "none")
#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_peripheral-prepost/predistancevs.preTIPS_bwsubj.pdf", plot=p,height=3, width=3.5)
ggsave("diversity_metrics_LCMS/beta_deicode_lcms_hepatic/predistancevs.preTIPS_bwsubj.pdf", plot=p,height=3, width=3.5)

pairwise.wilcox.test(perif_dist_pr$distance, perif_dist_pr$HE_pre_comb,
                     p.adjust.method="fdr")

#perif
# 0_0  
# 1_1 0.029

#hep
# 0_0 
# 1_1 0.29

#perif_dist_pr<-fread("diversity_metrics_LCMS/beta_deicode_lcms_peripheral-prepost/distance-matrix.tsv")%>%
perif_dist_pr<-fread("diversity_metrics_LCMS/beta_deicode_lcms_hepatic/distance-matrix.tsv")%>%
  gather(V2,distance,-V1) %>% 
  filter(V1 %in% preIds) %>%
  filter(V2 %in% preIds) %>%
  left_join(.,md1, by="V1") %>%
  left_join(.,md2, by="V2") %>%
  filter(!(V1_ATTRIBUTE_pt_id==V2_ATTRIBUTE_pt_id)) %>%
  mutate(HE_Wpost_comb=paste(V1_ATTRIBUTE_Worst_PostTIPS_HE_mod,V2_ATTRIBUTE_Worst_PostTIPS_HE_mod, sep="_")) %>%
  group_by(distance) %>% 
  filter(row_number()==1) %>%
  filter(HE_Wpost_comb %in% c("0_0","1_1","2+_2+"))

p<-ggplot(perif_dist_pr, aes(x=HE_Wpost_comb, y=distance,fill=HE_Wpost_comb)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#DC4D40","#754B27","#60B358"))+scale_fill_manual(values=c("#DC4D40","#754B27","#60B358"))+
  labs(x="HE grade",y="between subject robust aitchison", title="Pre TIPS")+
  theme(legend.position = "none")
#ggsave("diversity_metrics_LCMS/beta_deicode_lcms_peripheral-prepost/predistancevs.postTIPS_bwsubj.pdf", plot=p,height=3, width=3.5)
ggsave("diversity_metrics_LCMS/beta_deicode_lcms_hepatic/predistancevs.postTIPS_bwsubj.pdf", plot=p,height=3, width=3.5)

pairwise.wilcox.test(perif_dist_pr$distance, perif_dist_pr$HE_Wpost_comb,
                     p.adjust.method="fdr")

#perif
# 0_0  1_1 
# 1_1   0.12 -   
#   2+_2+ 0.12 0.98

#hep
# 0_0     1_1   
# 1_1   0.0018  -     
#   2+_2+ 4.1e-05 0.9542

##################################################
library(MASS)
library(ggord)
library(klaR)

md<-fread("/mnt/zarrinpar/scratch/sfloresr/HE/HE_metab/data_files/metadata_table/metadata_table-Perif-PrePost.txt") %>%
  select(sampleid,ATTRIBUTE_blood_procedure)
mtb_counts<-fread("data_files/LCMS/quantification_table/quantification_table-peripheral-prepost.txt") %>%
  mutate(`row ID`=paste("s",`row ID`,sep="")) %>%
  column_to_rownames("row ID") %>% t() %>% as.data.frame() %>% rownames_to_column("sampleid") %>%
  left_join(.,md,by="sampleid") %>%column_to_rownames("sampleid")%>%select(-c(429,458,557,558,586,588,590,595))

#mtb_counts<-mtb_counts[, colSums(mtb_counts != 0) > 0]

set.seed(12)
ind <- sample(2, nrow(mtb_counts),
              replace = TRUE,
              prob = c(0.6, 0.4))
training <- mtb_counts[ind==1,]
testing <- mtb_counts[ind==2,]

linear <- lda(ATTRIBUTE_blood_procedure~., training)
linear

p <- predict(linear, training)
ldahist(data = p$x[,1], g = training$ATTRIBUTE_blood_procedure)#

#ggord(linear, training$ATTRIBUTE_blood_procedure) only have LD1

#partimat(ATTRIBUTE_blood_procedure~., data = training, method = "lda")
