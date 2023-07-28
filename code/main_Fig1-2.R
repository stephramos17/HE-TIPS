#Main Code for the analysis of Figure 1-2

#Libraries
library("tidyverse")
library("data.table")
library("qiime2R")
library("ggpubr")
library("ggsci")
library("phyloseq")
library("vegan")

# Set directory
setwd("~/Notebooks/sfloresr/HE/HE-TIPS")

# subset counts table and metadata -------------------

# subsetting hepatic vein data
mdH<-fread("data/gnps_output/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-8068b7cb-view_all_clusters_withID/metadata_table/metadata_table-00000.tsv")%>%rename(sample_name=filename)%>%
  filter(ATTRIBUTE_blood_origin=="hepatic")%>%rename(sampleid=ATTRIBUTE_sampleID)%>%
  select(sampleid,everything())%>%select(-sample_name)
mdH$ATTRIBUTE_blood_procedure<-factor(mdH$ATTRIBUTE_blood_procedure,c("pre","post","bd","return"))

write.table(mdH,"data/processed/metadata_table-hepatic.txt",sep = "\t",row.names = FALSE,quote=FALSE)

HEdf_hepatic<-fread("data/gnps_output/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-8068b7cb-view_all_clusters_withID/quantification_table/quantification_table-00000.csv")%>%
  select(`row ID`, all_of(mdH$sampleid))
write.table(HEdf_hepatic,"data/processed/quantification_table-hepatic.txt",sep = "\t",row.names = FALSE,quote=FALSE)

# subsetting peripheral vein data

mdP<-fread("data/gnps_output/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-8068b7cb-view_all_clusters_withID/metadata_table/metadata_table-00000.tsv")%>%dplyr::rename(sample_name=filename)%>%
  filter(ATTRIBUTE_blood_origin=="peripheral",ATTRIBUTE_blood_procedure!="het")%>%
  dplyr::rename(sampleid=ATTRIBUTE_sampleID)%>%
  select(sampleid,everything())%>%select(-sample_name)
mdP$ATTRIBUTE_blood_procedure<-factor(mdP$ATTRIBUTE_blood_procedure,c("pre","post","bd","return"))

write.table(mdP,"data/processed//metadata_table-peripheral.txt",sep = "\t",row.names = FALSE,quote=FALSE)

HEdf_Perif<-fread("data/gnps_output/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-8068b7cb-view_all_clusters_withID/quantification_table/quantification_table-00000.csv")%>%
  select(`row ID`, all_of(mdP$sampleid))
write.table(HEdf_Perif,"data/processed/quantification_table-peripheral.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#subsetting post hepatic vein data

mdPoh<-fread("data/gnps_output/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-8068b7cb-view_all_clusters_withID/metadata_table/metadata_table-00000.tsv")%>%dplyr::rename(sample_name=filename)%>%
  filter(ATTRIBUTE_blood_procedure=="post" & ATTRIBUTE_blood_origin=="hepatic")%>%
  dplyr::rename(sampleid=ATTRIBUTE_sampleID)%>%select(sampleid,everything())%>%select(-sample_name)%>%
  mutate(WPostTIPS_HE=ifelse(ATTRIBUTE_Worst_PostTIPS_HE==3|ATTRIBUTE_Worst_PostTIPS_HE==2,"2+",ATTRIBUTE_Worst_PostTIPS_HE))
mdPoh$ATTRIBUTE_blood_procedure<-factor(mdPoh$ATTRIBUTE_blood_procedure,c("pre","post","bd","return"))
mdPoh$ATTRIBUTE_blood_origin<-factor(mdPoh$ATTRIBUTE_blood_origin,c("peripheral","hepatic"))

write.table(mdPoh,"data/processed/metadata_table-post-hepatic.txt",sep = "\t",row.names = FALSE,quote=FALSE)

HEdf_PostH<-fread("data/gnps_output/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-8068b7cb-view_all_clusters_withID/quantification_table/quantification_table-00000.csv")%>%
  select(`row ID`, all_of(mdPoh$sampleid))
write.table(HEdf_PostH,"data/processed/quantification_table-post-hepatic.txt",sep = "\t",row.names = FALSE,quote=FALSE)

# subsetting peripheral vein data --just pre vs. post

mdPP<-fread("data/gnps_output/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-8068b7cb-view_all_clusters_withID/metadata_table/metadata_table-00000.tsv")%>%dplyr::rename(sample_name=filename)%>%
  filter(ATTRIBUTE_blood_origin=="peripheral",ATTRIBUTE_blood_procedure!="het")%>%
  dplyr::rename(sampleid=ATTRIBUTE_sampleID)%>%
  select(sampleid,everything())%>%select(-sample_name)
mdPP$ATTRIBUTE_blood_procedure<-factor(mdPP$ATTRIBUTE_blood_procedure,c("pre","post","bd","return"))

write.table(mdPP,"data/processed/metadata_table-peripheral-prepost.txt",sep = "\t",row.names = FALSE,quote=FALSE)

HEdf_PerifPP<-fread("data/gnps_output/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-8068b7cb-view_all_clusters_withID/quantification_table/quantification_table-00000.csv")%>%
  select(`row ID`, all_of(mdPP$sampleid))
write.table(HEdf_PerifPP,"data/processed/quantification_table-peripheral-prepost.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#cleaning mtb file -------------------
mtbid_dict<-fread("data/gnps_output/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-8068b7cb-view_all_clusters_withID/DB_result/32405003deff48dfb9e7cb571ecefc23.tsv")%>%
  rename(rowID=`#Scan#`)

md<-fread("data/gnps_output/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-8068b7cb-view_all_clusters_withID/metadata_table/metadata_table-00000.tsv")%>%rename(sample_name=filename)
md$ATTRIBUTE_blood_procedure<-factor(md$ATTRIBUTE_blood_procedure,c("pre","post","bd","return"))

mtb_all<-fread("data/gnps_output/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-8068b7cb-view_all_clusters_withID/quantification_table/quantification_table-00000.csv")%>%
  gather(ATTRIBUTE_sampleID,counts, -`row ID`,-`row m/z`,-`row retention time`)%>% rename(rowID=`row ID`)%>%
  mutate(log_counts=log10(counts+1))%>% left_join(.,md, by="ATTRIBUTE_sampleID")%>%
  left_join(.,mtbid_dict,by="rowID")

write.table(mtb_all,"data/processed/HEmtbdata_allwmd.txt",sep = "\t",row.names = FALSE,quote=FALSE)

# paired wilcoxon analysis -------------------
##peripheral vein
mtb_all<-fread("data/processed/HEmtbdata_allwmd.txt")

perifPrPo<-run_paired_wilcox(mtb_all,"pre_vs_post","peripheral")
perifPrbd<-run_paired_wilcox(mtb_all,"pre_vs_bd","peripheral")
perifPobd<-run_paired_wilcox(mtb_all,"post_vs_bd","peripheral")

perif_pvals<-rbind(perifPrPo,perifPrbd,perifPobd)
write.table(perif_pvals,"data/results/Perif_HEVisitvslogMtbAbun_PAIREDwilcoxFDR_allpvals.txt",sep = "\t",
            row.names = FALSE,quote=FALSE)
perif_sigpvals<-perif_pvals%>%filter(padj<0.1)
write.table(perif_sigpvals,"data/results/Perif_HEVisitvslogMtbAbun_PAIREDwilcoxFDR0.1_sigpvals.txt",sep = "\t",
            row.names = FALSE,quote=FALSE)

##hepatic vein
hepPrPo<-run_paired_wilcox(mtb_all,"pre_vs_post","hepatic")
write.table(hepPrPo,"data/results/hepatic_HEVisitvslogMtbAbun_PAIREDwilcoxFDR_allpvals.txt",sep = "\t",
            row.names = FALSE,quote=FALSE)
hepPrPo_sigpvals<-hepPrPo%>%filter(padj<0.1)
write.table(hepPrPo_sigpvals,"data/results/hepatic_HEVisitvslogMtbAbun_PAIREDwilcoxFDR0.1_sigpvals.txt",sep = "\t",
            row.names = FALSE,quote=FALSE)

# summarise paired wilcoxon analysis -------------------
##heatmap for peripheral and hepatic (pre vs post)
sigmtbsP<-fread("data/results/Perif_HEVisitvslogMtbAbun_PAIREDwilcoxFDR0.1_sigpvals.txt")%>%
  mutate(ATTRIBUTE_blood_origin="peripheral")%>%filter(comparison=="pre_vs_post")
sigmtbsH<-fread("data/results/hepatic_HEVisitvslogMtbAbun_PAIREDwilcoxFDR0.1_sigpvals.txt")%>%
  mutate(ATTRIBUTE_blood_origin="hepatic")
sigmtbs<-rbind(sigmtbsP,sigmtbsH)%>%rename(rowID=mtb_id)
write.table(sigmtbs,"data/results/ALLbldorg_HEVisitvslogMtbAbun_PAIREDsigwilcoxFDR0.1_prepost.txt",sep = "\t",
            row.names = FALSE,quote=FALSE)

sigmtbs_annotA<-sigmtbs%>%rename(ATTRIBUTE_blood_procedure=comparison)%>%
  mutate(rowID=as.numeric(rowID),ATTRIBUTE_blood_procedure=ifelse(ATTRIBUTE_blood_procedure=="pre_vs_post","post",ATTRIBUTE_blood_procedure))%>%
  mutate(ATTRIBUTE_blood_procedure=ifelse(ATTRIBUTE_blood_procedure=="pre_vs_bd","bd",ATTRIBUTE_blood_procedure))%>%
  filter(ATTRIBUTE_blood_procedure=="post"|ATTRIBUTE_blood_procedure=="bd")%>%
  mutate(stars=cut(padj, breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1,Inf), label=c("***", "**", "*", ".","")))

mtb_allA<-fread("data/processed/HEmtbdata_allwmd.txt")%>%
  filter(ATTRIBUTE_blood_origin=="peripheral"|ATTRIBUTE_blood_origin=="hepatic")%>%
  filter(ATTRIBUTE_blood_procedure=="pre"|ATTRIBUTE_blood_procedure=="post")

sigmtb_cleanA<-mtb_allA%>%filter(rowID %in% sigmtbs$rowID)%>% left_join(.,sigmtbs_annotA, by=c("rowID","ATTRIBUTE_blood_procedure","ATTRIBUTE_blood_origin"))%>%
  mutate(label_name=ifelse(!is.na(Compound_Name),paste("mtb ",rowID," ",Compound_Name,sep=""),paste("mtb ",rowID," m/z",`row m/z`)),
         rowID=as.numeric(rowID))%>%group_by(label_name,ATTRIBUTE_blood_procedure,ATTRIBUTE_blood_origin,stars)%>%summarise(mn_log_counts=mean(log_counts))%>%
  mutate(ATTRIBUTE_blood_origin=factor(ATTRIBUTE_blood_origin, levels = c("peripheral", "hepatic")),
         ATTRIBUTE_blood_procedure=factor(ATTRIBUTE_blood_procedure, levels = c("pre", "post","bd")))

write.table(sigmtb_cleanA,"data/results/prevpost_PerifHep_heatmapdata.txt",sep = "\t",
            row.names = FALSE,quote=FALSE)

##heatmap for peripheral (pre vs post vs bd)

sigmtbsPh<-fread("data/results/Perif_HEVisitvslogMtbAbun_PAIREDwilcoxFDR0.1_sigpvals.txt")%>%
  mutate(ATTRIBUTE_blood_origin="peripheral")%>%rename(rowID=mtb_id)

sigmtbs_annotB<-sigmtbsPh%>%rename(ATTRIBUTE_blood_procedure=comparison)%>%
  mutate(rowID=as.numeric(rowID),ATTRIBUTE_blood_procedure=ifelse(ATTRIBUTE_blood_procedure=="pre_vs_post","post",ATTRIBUTE_blood_procedure))%>%
  mutate(ATTRIBUTE_blood_procedure=ifelse(ATTRIBUTE_blood_procedure=="pre_vs_bd","bd",ATTRIBUTE_blood_procedure))%>%
  filter(ATTRIBUTE_blood_procedure=="post"|ATTRIBUTE_blood_procedure=="bd")%>%
  mutate(stars=cut(padj, breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ".","")))

mtb_allB<-fread("data/processed/HEmtbdata_allwmd.txt")%>%
  filter(ATTRIBUTE_blood_origin=="peripheral")%>%
  filter(ATTRIBUTE_blood_procedure=="pre"|ATTRIBUTE_blood_procedure=="post"|ATTRIBUTE_blood_procedure=="bd")

sigmtb_cleanB<-mtb_allB%>%filter(rowID %in% sigmtbsP$rowID)%>% left_join(.,sigmtbs_annotB, by=c("rowID","ATTRIBUTE_blood_procedure","ATTRIBUTE_blood_origin"))%>%
  mutate(label_name=ifelse(!is.na(Compound_Name),paste("mtb ",rowID," ",Compound_Name,sep=""),paste("mtb ",rowID," m/z",`row m/z`)),
         rowID=as.numeric(rowID))%>%filter(!is.na(Compound_Name))%>%
  group_by(label_name,ATTRIBUTE_blood_procedure,ATTRIBUTE_blood_origin,stars)%>%summarise(mn_log_counts=mean(log_counts))%>%
  mutate(ATTRIBUTE_blood_origin=factor(ATTRIBUTE_blood_origin, levels = c("peripheral", "hepatic")),
         ATTRIBUTE_blood_procedure=factor(ATTRIBUTE_blood_procedure, levels = c("pre", "post","bd")))

write.table(sigmtb_cleanB,"data/results/prevpostvsbd_PerifOnly_heatmapdata.txt",sep = "\t",
            row.names = FALSE,quote=FALSE)
####################################################################
#                  Figure 1 - mtb changes over TIPS procedure          
####################################################################

# 1.B. RPCA plot of hepatic vein samples--------------
hep_ord <- read_qza("data/results/ordination-hepatic.qza") #check deicode.sh

rpcaH<-hep_ord$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  rename(sampleid=SampleID)%>%
  left_join(mdH,by="sampleid")

pH<-rpcaH %>%
  ggplot(aes(x=PC1, y=PC2, color=ATTRIBUTE_blood_procedure)) +
  geom_point(alpha=1.0, shape=15) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +stat_ellipse(type = "t", linetype = 2,aes(group = ATTRIBUTE_blood_procedure))+
  scale_color_manual(values=c("#52A7DC","#DF643B","#2D653D"))+
  labs(x =paste("PC1 (",round(hep_ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(hep_ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       color="TIPS visit")+ggtitle("Hepatic")+ theme(plot.title = element_text(face = "bold"))

ggsave("data/figures/PCoAvisitHVOnlyAnnot_wellipse.pdf", plot=pH,height=3, width=3)

# 1.C. RPCA plot of peripheral vein samples--------------
perif_ord <- read_qza("data/results/ordination-peripheral.qza") #check deicode.sh

rpcaP<-perif_ord$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  dplyr::rename(sampleid=SampleID)%>%
  left_join(mdP,by="sampleid")

pP<-rpcaP %>%
  ggplot(aes(x=PC1, y=PC2, color=ATTRIBUTE_blood_procedure)) +
  geom_point(alpha=1.0) + 
  theme_pubr() +stat_ellipse(type = "t", linetype = 2,aes(group = ATTRIBUTE_blood_procedure))+
  scale_color_manual(values=c("#52A7DC","#DF643B","#2D653D"))+
  labs(x =paste("PC1 (",round(perif_ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(perif_ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       color="TIPS visit")+ggtitle("Peripheral")+ theme(plot.title = element_text(face = "bold"))

ggsave("data/figures/PCoAvisitPerifOnlyAnnot_wellipse.pdf", plot=pP,height=3, width=3.1)

# 1.D. Natural log ratio of bile acids over the TIPS procedure in peripheral vein--------------
df_div<-fread("data/results/sample_plot_data_bab_flip.tsv")%>%select(1:3)%>%
  mutate(ATTRIBUTE_blood_procedure=factor(ATTRIBUTE_blood_procedure, levels = c("pre", "post","bd")))

Natlog_p<-ggplot(df_div, aes(x=ATTRIBUTE_blood_procedure, y=Current_Natural_Log_Ratio,fill=ATTRIBUTE_blood_procedure)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#52A7DC","#DF643B","#2D653D"))+
  scale_fill_manual(values=c("#52A7DC","#DF643B","#2D653D"))+
  labs(x="TIPS visit",y="Natural Log Ratio")+
  theme(legend.position = "none")

ggsave("data/figures/deicode_bab_natlog.pdf",Natlog_p, height=3, width=3.5)

pairwise.wilcox.test(df_div$Current_Natural_Log_Ratio, df_div$ATTRIBUTE_blood_procedure,
                     p.adjust.method="fdr")
 
#         pre    post  
# post  0.0020    -     
#  bd   0.7180  0.0032

# 1.E. heatmap of significantly expressed mtb hits--------------

# 1.F. example: differential expression of GUDCA pre vs. post TIPS (peripheral) --------------
mtb_all<-fread("data/processed/HEmtbdata_allwmd.txt")
pval_dfF<-fread("data/results/Perif_HEVisitvslogMtbAbun_PAIREDwilcoxFDR0.1_sigpvals.txt")%>%
  filter(comparison=="pre_vs_post")
pltF<-paired_plots(mtb_all,c(515),"pre_vs_post","peripheral",pval_dfF,1,c("#52A7DC","#DF643B"))
ggsave("data/figures/mtb382515_Perifprepost_Paired.pdf", plot=pltF)

# 1.G. example: differential expression of bile acids post vs. bd TIPS (peripheral) --------------
mtb_all<-fread("data/processed/HEmtbdata_allwmd.txt")
pval_dfG<-fread("data/results/Perif_HEVisitvslogMtbAbun_PAIREDwilcoxFDR0.1_sigpvals.txt")%>%
  filter(comparison=="post_vs_bd")
pltG<-paired_plots(mtb_all,c(176,177,181),"post_vs_bd","peripheral",pval_dfG,3,c("#DF643B","#2D653D"))
ggsave("diff_HE_analysis/BA_Perifpostbd_Paired_shortannot.pdf", plot=pltG)

####################################################################
#                  Figure 2 - mtb changes associated with HE grade        
####################################################################

# 2.A. RPCA of post hepatic vein samples showing HE grade--------------
grd_ord <- read_qza("data/results/ordination-post-hepatic.qza")

rpcaPoH<-grd_ord$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  dplyr::rename(sampleid=SampleID)%>%
  left_join(mdH,by="sampleid")

poH<-rpcaPoH %>%
  ggplot(aes(x=PC1, y=PC2, color=ATTRIBUTE_Worst_PostTIPS_HE_mod,shape=ATTRIBUTE_blood_origin)) +
  geom_point(alpha=1.0) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +stat_ellipse(type = "t", linetype = 2,aes(group = ATTRIBUTE_Worst_PostTIPS_HE_mod))+
  scale_shape_manual(values=c(15,16,17),name="Blood Origin") +
  scale_color_manual(values=c("#39B54A", "#283891", "#EF3E36"))+
  labs(x =paste("PC1 (",round(grd_ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(grd_ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       color="HE grade")+ggtitle("PostTIPS") + theme(plot.title = element_text(face = "bold"), legend.position = "left")
ggsave("data/figures/PCoAvisitPostOnlyHEgradeHep_wellipse.pdf", plot=poH,height=3, width=4.5)

##0
rpca_0<-rpcaPoH%>%filter(ATTRIBUTE_Worst_PostTIPS_HE_mod=="0")%>%dplyr::select(2:3)
eig_dat <- eigen(cov(rpca_0))$values 
vec <- sqrt(5.991* eig_dat)
pi * vec[1] * vec[2] 
#A=0.0577387

##1
rpca_1<-rpcaPoH%>%filter(ATTRIBUTE_Worst_PostTIPS_HE_mod=="1")%>%dplyr::select(2:3)
eig_dat <- eigen(cov(rpca_1))$values 
vec <- sqrt(5.991* eig_dat)
pi * vec[1] * vec[2] 
#A=1.128352

##2+
rpca_2p<-rpcaPoH%>%filter(ATTRIBUTE_Worst_PostTIPS_HE_mod=="2+")%>%dplyr::select(2:3)
eig_dat <- eigen(cov(rpca_2p))$values 
vec <- sqrt(5.991* eig_dat)
pi * vec[1] * vec[2] 
#A=1.191518

# 2.B. change in portal pressure over TIPS study --------------


calculate_portal_gradient <- function(ptFile, pgFile) {
  
  # The goal of this function is to plot gradient change (Fig 2b).
  
  # We do so by using the direct portal vein pressure in mmHg 
  # before and after TIPS.
  
  # Parameters:
  # ----------
  # ptFile: patient metadata file path
  # pgFile: patient gradient pressure 
  
  # Returns
  # -------
  
  
  # Read input files and merde by id
  ptdat <- readxl::read_excel(ptFile, na = c("", " ", "NA"))
  pgdat <- read.csv(pgFile)
  dat <- merge(ptdat, pgdat, by="pt_id")
  
  # Calculate gradient
  grad <- dat %>%
    subset(select=c("Worst_PostTIPS_HE_mod",  "Pre.PV.gradient", "Post.PV.gradient")) %>%
    drop_na() %>%
    mutate(`pressure_change` = Pre.PV.gradient - Post.PV.gradient )
  
  # Set colors for plot
  he_colors <- c("#39B54A", "#283891" , "#EF3E36")
  
  # Plot data
  p <- ggplot(grad, aes(x=Worst_PostTIPS_HE_mod, y=pressure_change, fill=Worst_PostTIPS_HE_mod)) +
    theme_bw()+
    geom_boxplot(outlier.shape = NA, alpha=0.4) +
    geom_dotplot(binaxis='y', stackdir='center', alpha=0.7,
                 position=position_dodge(1)) +
    xlab("HE grading") +
    ylab("mmHg change [mmHgpre - mmHgpost]")+
    theme(legend.position="bottom") +
    scale_fill_manual(values=he_colors)+
    stat_compare_means(method = "kruskal", label = "p.format",
                       label.x.npc = 0.3, label.y.npc = 0.9)
  ggsave("./figures/gradient_change.pdf", width=4, height=7)
  ggsave("./figures/gradient_change.png", width=4, height=7)
  
  
}

calculate_portal_gradient(ptFile="./data/processed/patient_overview_short.xlsx",
                          pgFile="./data/processed/TIPS_pressure_data_20230310.csv")


# 2.C. within individual changes of metabolite dissimilarity in peripheral vein--------------
mdP<-fread("data/processed/metadata_table-peripheral-prepost.txt") %>%
  dplyr::select(sampleid,ATTRIBUTE_pt_id, ATTRIBUTE_blood_procedure,
                ATTRIBUTE_PreTIPS_HE_mod,ATTRIBUTE_Worst_PostTIPS_HE_mod)%>%
  spread(ATTRIBUTE_blood_procedure,sampleid)

perif_dist<-fread("data/results/beta_deicode_lcms_peripheral-prepost/distance-matrix.tsv")%>%
  gather(V2,distance,-V1) %>%dplyr::rename(pre=V1, post=V2) %>%
  right_join(.,mdP,by=c("pre","post")) %>%
  drop_na() %>%
  mutate(ATTRIBUTE_PreTIPS_HE_mod=as.factor(ATTRIBUTE_PreTIPS_HE_mod),
         ATTRIBUTE_Worst_PostTIPS_HE_mod=as.factor(ATTRIBUTE_Worst_PostTIPS_HE_mod))
#ranked plot
pwP<-ggplot(perif_dist, aes(x=fct_reorder(ATTRIBUTE_pt_id,distance), y=distance, fill=ATTRIBUTE_Worst_PostTIPS_HE_mod)) +
  geom_bar(stat="identity")+theme_classic() + scale_color_manual(values=c("#39B54A", "#283891", "#EF3E36"))+
  scale_fill_manual(values=c("#39B54A", "#283891", "#EF3E36"))+
  labs(x="subject ID",y="robust aitchison distance", title="Post TIPS grading", fill="HE grading")+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("data/figures/peripheral_distancepost_bysubject.pdf", plot=pwP,height=3.5, width=5)

#boxplot
pwPs<-ggplot(perif_dist, aes(x=ATTRIBUTE_Worst_PostTIPS_HE_mod, y=distance,fill=ATTRIBUTE_Worst_PostTIPS_HE_mod)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#39B54A", "#283891", "#EF3E36"))+
  scale_fill_manual(values=c("#39B54A", "#283891", "#EF3E36"))+
  labs(x="HE grade",y="robust aitchison distance", title="Post TIPS grading")+
  theme(legend.position = "none")
ggsave("data/figures/peripheral_distancevs.postTIPS.pdf", plot=pwPs,height=3, width=3.5)

pairwise.wilcox.test(perif_dist$distance, perif_dist$ATTRIBUTE_Worst_PostTIPS_HE_mod,
                     p.adjust.method="fdr")

#        0    1   
#    1  0.87  -   
#   2+ 0.87 0.87

# 2.D. within individual changes of metabolite dissimilarity in hepatic vein--------------
mdH<-fread("data/processed/metadata_table-hepatic.txt") %>%
  dplyr::select(sampleid,ATTRIBUTE_pt_id, ATTRIBUTE_blood_procedure,
                ATTRIBUTE_PreTIPS_HE_mod,ATTRIBUTE_Worst_PostTIPS_HE_mod)%>%
  spread(ATTRIBUTE_blood_procedure,sampleid)

hepatic_dist<-fread("data/results/beta_deicode_lcms_hepatic/distance-matrix.tsv")%>%
  gather(V2,distance,-V1) %>%dplyr::rename(pre=V1, post=V2) %>%
  right_join(.,mdH,by=c("pre","post")) %>%
  drop_na() %>%
  mutate(ATTRIBUTE_PreTIPS_HE_mod=as.factor(ATTRIBUTE_PreTIPS_HE_mod),
         ATTRIBUTE_Worst_PostTIPS_HE_mod=as.factor(ATTRIBUTE_Worst_PostTIPS_HE_mod))
#ranked plot
pwH<-ggplot(hepatic_dist, aes(x=fct_reorder(ATTRIBUTE_pt_id,distance), y=distance, fill=ATTRIBUTE_Worst_PostTIPS_HE_mod)) +
  geom_bar(stat="identity")+theme_classic() + scale_color_manual(values=c("#39B54A", "#283891", "#EF3E36"))+
  scale_fill_manual(values=c("#39B54A", "#283891", "#EF3E36"))+
  labs(x="subject ID",y="robust aitchison distance", title="Post TIPS grading", fill="HE grading")+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("data/figures/hepatic_distancepost_bysubject.pdf", plot=pwH,height=3.5, width=5)

#boxplot
pwHs<-ggplot(hepatic_dist, aes(x=ATTRIBUTE_Worst_PostTIPS_HE_mod, y=distance,fill=ATTRIBUTE_Worst_PostTIPS_HE_mod)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#39B54A", "#283891", "#EF3E36"))+
  scale_fill_manual(values=c("#39B54A", "#283891", "#EF3E36"))+
  labs(x="HE grade",y="robust aitchison distance", title="Post TIPS grading")+
  theme(legend.position = "none")
ggsave("data/figures/hepatic_distancevs.postTIPS.pdf", plot=pwHs,height=3, width=3.5)

pairwise.wilcox.test(hepatic_dist$distance, hepatic_dist$ATTRIBUTE_Worst_PostTIPS_HE_mod,
                     p.adjust.method="fdr")

#        0     1    
#    1  0.027  -    
#   2+ 0.036 0.272
