# main_Fig3-5.R
# Author: Carolina Dantas
# Date: Mar 9 2023
# Purpose: code to perform analysis of metabolomics data and generate Figures 3-5

# Libraries
library(tidyverse)
library(ggpubr)
library(viridis)
library(rstatix)

# Set directories and analysis paths
#setwd("~/Dropbox/ucsd/projects/he/analysis/HE-TIPS/")
fbmnDir <- "./data/gnps_output/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-8068b7cb-view_all_clusters_withID/"
analysis_name <-"20230309_"
source("./code/function_lcms.R")

# Set dir and paths
inDir <- "./data/gnps_output" 
file_mdat <- file.path(inDir, "filtered_results_20220430.csv" )
chemdir_files <- list.files(inDir, pattern="CHEMDIR.*tsv", recursive=TRUE, full.names = T)
fig_pre <- paste0("./figures/", analysis_name)
res_pre <- paste0("./results/", analysis_name)

# Use selected colors
my_colors <- list(pre="#52A7DC", post="#DF643B", hep="#EE207C", per="#2E3192")
he_colors <- c("#39B54A", "#283891" , "#EF3E36")


# Load and process input data -------------------

# Read input data
qtab <- load_quant_tab_lcms() # quantification table
db <- load_db_tab_lcms() # main hits data table
cdat <- load_clusterinfo_tab_lcms() # clusterinfo data table
mdat <- load_metadata_lcms() # metadata table
mdb <- add_cluster_subclass_db() # get db with component index and subclass


# Pre-process data
lqtab <- get_long_qtab_w_mdat()
lqtab <- add_log10_lqtab(dat=lqtab)
lqtab

# Get paired hepatic data heppair_lqtab
heppair_lqtab <- subset(lqtab, Sample %in%  c("Pre-HV", "Post-HV")) %>%
  get_paired_dat_pt()

# Get paired  peripheral data
perpair_lqtab <- subset(lqtab, Sample %in%  c("Pre-P", "Post-P")) %>%
  get_paired_dat_pt()

####################################################################
#                  Figure 3 - Change from baseline          
####################################################################


# Pre-process data for change from baseline analysis --------------

# Get change from baseline for peripheral data
per_change <- get_change_from_baseline(dat=perpair_lqtab) %>%
  subset(., select = c("row.ID", "componentindex", "log10_dpeak", "Worst_PostTIPS_HE_mod", "blood_origin"))


# Get change form baseline for hepatic data
hep_change <- get_change_from_baseline(dat=heppair_lqtab) %>%
  subset(., select = c("row.ID", "componentindex", "log10_dpeak", "Worst_PostTIPS_HE_mod", "blood_origin"))

# Combined
all_change <- rbind(per_change, hep_change)


write.csv(all_change, file = paste0(res_pre, "fig3-table.csv"), quote=F, row.names = F)

# Save tables

# Plot Change from baseline results

# 3.B. Density plot for change from baseline ------------
dir.create("./figures")
dir.create("./results")

# Peripheral only ---------

make_dens_plot(per_change)
ggsave(paste0(fig_pre, "dens_per_deltapeak.pdf"))
run_ks_for_groups(per_change)
# 0 vs 1 P=0.0009974
# 0 vs 2+ P = 1.255e-13
# 1 vs 2+ P = 2.429e-09


# Hepatic only -------

make_dens_plot(hep_change)
ggsave(paste0(fig_pre, "dens_hep_deltapeak.pdf"))
run_ks_for_groups(hep_change)
# 0 vs 1 P=0.08188
# 0 vs 2+ P = 5.71e-08
# 1 vs 2+ P = < 2.2e-16


# 3.C. Main significant clusters with delta changes --------------


# Identify significant clusters ------


# Peripheral --------
get_significance_clusters(per_change) %>%
  subset(p.adj < 0.05)
# c(1, 8, 12)

# Hepatic -------
get_significance_clusters(hep_change) %>%
  subset(p.adj < 0.05)
# c(5, 8, 19)


get_significance_clusters(per_change) %>%
  add_significance("p.adj") %>%
  write.csv(file="./results/kruskal_log10dpeak_HEgrading_per_sig_clusters.csv", quote=F, row.names = F)

get_significance_clusters(hep_change) %>%
  add_significance("p.adj") %>%
  write.csv(file="./results/kruskal_log10dpeak_HEgrading_heep_sig_clusters.csv", quote=F, row.names = F)


# Violin plot of change from baseline -------

# Peripheral
per_change %>%
  subset(componentindex %in% c(1, 8, 12)) %>%
  make_violinplot_byGrade()
ggsave(paste0(fig_pre, "dviolin_clusters_deltapeak_per.pdf"), width = 9, height = 3.5)

per_change %>%
  subset(componentindex %in% c(1, 8, 12)) %>%
  compare_means(formula = log10_dpeak ~ Worst_PostTIPS_HE_mod , method = "wilcox", 
                group.by = "componentindex",
                comparisons = list(c("0", "1"), c("0", "2+"), c("1", "2+"))) %>%
  add_significance("p.adj") %>%
  write.csv(file="./results/wilcox_log10dpeak_HEgrading_per_sig_clusters.csv", 
            quote=F, row.names = F)

# Hepatic
hep_change %>%
  subset(componentindex %in% c(5, 8, 19)) %>%
  make_violinplot_byGrade()
ggsave(paste0(fig_pre, "dviolin_clusters_deltapeak_hep.pdf"), width = 9, height = 3.5)


hep_change %>%
  subset(componentindex %in% c(5, 8, 19)) %>%
  compare_means(formula = log10_dpeak ~ Worst_PostTIPS_HE_mod , method = "wilcox", 
                group.by = "componentindex",
                comparisons = list(c("0", "1"), c("0", "2+"), c("1", "2+"))) %>%
  add_significance("p.adj") %>%
  write.csv(file="./results/wilcox_log10dpeak_HEgrading_hep_sig_clusters.csv", 
            quote=F, row.names = F)


# Plot main clusters
all_change %>%
  subset(componentindex %in% c("8")) %>%
  make_violinplot_byGrade() +
  facet_wrap(~blood_origin) -> main_bileacid_panel
main_bileacid_panel
ggsave(paste0(fig_pre, "dviolin_clusters_deltapeak_bileacid_main.pdf"))


per_change %>%
  subset(componentindex == "1") %>%
  make_violinplot_byGrade() -> main_gpc_panel


# 3.A. Barplot with ordered metabolite--------------

# Get sorted mean values for each metabolite
sm_per <- get_sorted_mb_4plot(dat=per_change, grade=0)
sm_hep <- get_sorted_mb_4plot(dat=hep_change, grade=0)

sm_all <- get_sorted_mb_4plot(transform(
  rbind(hep_change, per_change), blood_origin="all"), grade=0)

# Plot segment by order in 0
sm_per %>%
  ggplot(aes(y=mean_change, x=row.ID)) +
  theme_pubr()+
  geom_segment( aes(x=row.ID, xend=row.ID, y=0, 
                    yend=mean_change, 
                    color=Worst_PostTIPS_HE_mod), alpha=0.7) +
  scale_color_manual(values=c("#39B54A", "#283891", "#EF3E36"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom") +
  xlab("Feature ranking")+
  facet_wrap(~blood_origin )-> segm_plot1
segm_plot1
ggsave(filename = paste0(fig_pre, "fig3_per_baselinechange_byHEgrading.pdf"),
       width=7, height = 3)

sm_hep %>%
  ggplot(aes(y=mean_change, x=row.ID)) +
  theme_pubr()+
  geom_segment( aes(x=row.ID, xend=row.ID, y=0, 
                    yend=mean_change, 
                    color=Worst_PostTIPS_HE_mod), alpha=0.7) +
  scale_color_manual(values=c("#39B54A", "#283891", "#EF3E36"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom") +
  xlab("Feature ranking")+
  facet_wrap(~blood_origin ) -> segm_plot1b
segm_plot1b
ggsave(filename = paste0(fig_pre, "fig3_hep_baselinechange_byHEgrading.pdf"),
       width=7, height = 3)


# Save  main panels together
panel_clusters  
ggarrange(plotlist = list(segm_plot1, segm_plot1b), 
          nrow = 2,  common.legend = T)
ggsave(filename = paste0(fig_pre, "_fig3_a_b.pdf"), width = 4, height=3 )

# Peripheral and Hepatic - facet
make_dens_plot(all_change) +
  facet_wrap(~blood_origin) +
  xlab("Change pre/post") +
  theme(panel.spacing = unit(2, "lines")) 
ggsave(paste0(fig_pre, "dens_per_hep_deltapeak.pdf"), width = 4.5, height=4  )


ggarrange(plotlist = list(main_bileacid_panel, main_gpc_panel), 
          nrow=1, widths = c(2,1.15), common.legend = T)
ggsave(paste0(fig_pre, "dviolin_clusters_deltapeak_main.pdf"), width = 7, height = 3)


####################################################################
#                 Figure 4 -  Bile Acids analysis
####################################################################



# Pre-process data ---------------

# Rename Sample for manuscript figures and 
# rename for consistency across participants 
lqtab <- rename_sample_abbreviations()

# Select all bile acids by cluster index
allbas <- subset(cdat, componentindex == 8, select=cluster.index)[[1]]
allbas

# Find significant BAs ----------------------

# No significance in Pre levels w/ kruskal

subset(lqtab, row.ID  %in% allbas & Sample == "Post-PIV" ) %>%
  run_stat_mbt_abundance(test="kruskal") %>%
  write.csv(file="./results/kruskal_abundance_bileacids_post-piv.csv", quote=F, row.names=F)

subset(lqtab, row.ID  %in% allbas & Sample == "Post-HV" ) %>%
  run_stat_mbt_abundance(test="kruskal") %>%
  write.csv(file="./results/kruskal_abundance_bileacids_post-hv.csv", quote=F, row.names=F)


# Significant  by wilcox b/e HE grades
subset(lqtab, row.ID  %in% allbas & Sample == "Post-PIV" ) %>%
  run_stat_mbt_abundance(test="wilcox") %>%
  subset(p <= 0.05)

subset(lqtab, row.ID  %in% allbas & Sample == "Post-HV" ) %>%
  run_stat_mbt_abundance(test="wilcox") %>%
  subset(p <= 0.05)


# Select significant bile acids for plots
# Most differences  are for post-PIV
subset(lqtab, row.ID  %in% allbas & Sample == "Post-PIV" ) %>%
  run_stat_mbt_abundance(test="kruskal") %>%
  subset(p <= 0.05, select=row.ID) %>% unlist() %>% as.vector() -> bas


# Plot specific BA levels ------------------
subset(lqtab, row.ID  %in% bas ) %>%
  subset(select=c("row.ID", "Sample", "pt_id", "log10", "Worst_PostTIPS_HE_mod")) %>%
  drop_na(Sample) %>%
  write.csv(., file = paste0(res_pre, "fig4-table.csv"), quote=F, row.names = F)


# Figure 4a - sig by HE grade post-PIV
subset(lqtab, row.ID  %in% bas & Sample == "Post-PIV" ) %>%
  plot_bileacids_violin() +
  facet_wrap(~row.ID) 
ggsave(paste0(fig_pre, "fig4a_bileacids.pdf"), width = 7, height = 2.5)


# Figure 4B-C - specific bile acids for each pt
subset(lqtab, row.ID  %in% bas ) %>%
  drop_na(Sample) %>%
  ggplot(aes(y=log10, x=Sample, group=1)) +
  theme_bw() +
  ylab("log10(abundance)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~row.ID) +
  stat_summary(fun = "mean", geom="line", color="grey", alpha=0)+
  stat_summary(data = ~ subset(.x, Worst_PostTIPS_HE_mod == "0"),
               fun = "mean", geom="line", color="black") +
  stat_summary(data = ~ subset(.x, Worst_PostTIPS_HE_mod == "0"),
               fun = "mean", geom="point", color="black") +
  stat_summary(data = ~ subset(.x, Worst_PostTIPS_HE_mod == "0"),
               fun.data = "mean_se", geom="ribbon", fill="black", alpha=0.3) -> lp_mean_bas

# Patient 11
lp_mean_bas +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_line(data = ~ subset(.x, pt_id %in% c("11_a")), color="red") +
  geom_point(data = ~ subset(.x, pt_id %in% c("11_a")), color="red")
ggsave(paste0(fig_pre, "fig4b_bileacids_pt11.pdf"),  width = 7, height = 3)

# Patient 12
lp_mean_bas +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_line(data = ~ subset(.x, pt_id %in% c("12_a")), color="orange") +
  geom_point(data = ~ subset(.x, pt_id %in% c("12_a")), color="orange")
ggsave(paste0(fig_pre, "fig4c_bileacids_pt12.pdf"),width = 7, height = 3)

####################################################################
#                 Figure 5 - chemprop analysis
####################################################################

# Read chemdir
cddat <- read_chemdir_files(chemdir_files)
selfloop <- load_selfloop_table(fbmnDir)

# Chemdir jobs
jobs <- read_chemdirjob_mdat(f=file_mdat)

# Add columns "clusters" and "abs_dmz"
df <- add_new_chemdir_cols(x=cddat)


# Merge selfloop and main chemdir output 
df <- merge(df, selfloop, all.x=TRUE,
            by=c("CLUSTERID1", "CLUSTERID2", "DeltaMZ", "Cosine", "ComponentIndex" ))

# Simplify naming
df$EdgeAnnotation <- as.factor(gsub("(C4H4N3O_cytosine).*", "\\1\\*", df$EdgeAnnotation))
df$EdgeAnnotation <- as.factor(gsub("(C6H10O5_C6H10O5:).*", "\\1\\*", df$EdgeAnnotation))


# Merge chemdir data and file metadata info
df <- merge(df, jobs[,c("group", "timeseries", "HE_status", "file_id", "condition",
                        "grading_value", "grading_timepoint")], 
            by="file_id") 

# Save combined results table
df %>%
  write.csv(., file = paste0(res_pre, "fig5-table.csv"), quote=F, row.names = F)

# Pre vs post - hepatic and peripheral -------------------

# Main delta m/z plots combined 

subset(df, !grepl("adduct", df$EdgeAnnotation) & 
         timeseries == "PreVsPost" &
         HE_status == "all" ) %>%
  plot_geompoint_dmz() +
  facet_wrap(~group, nrow=2) +
  xlim(c(-100, 100)) -> main_dmz
main_dmz
ggsave(paste0(fig_pre,  "pointplot_deltaMZ_all_pubr.pdf"), width = 7, height = 4)

cluster_ids <- c("33_190", "219_733", "527_733", "29_406")
subset(df, timeseries == "PreVsPost" &
         HE_status == "all" & clusters %in% cluster_ids )



# Do analysis based on worst HE grade post TIPS ------------------------------

# Select chemdir > 2 by cluster numbers and plot those -------------
grad <- subset(df, grading_timepoint =="wpostTIPS") 
grad

# Select clusters
selec_clusters <- unique(df[with(df, max_chemdir > 2) ,]$clusters)


# Barplot of chemprop main differences -------------------


#  Figure  5e - barplot
# Barplot chemprop main parirs # 

subset(grad, group %in% c("hep")) %>%
  get_diff_chemdir() %>%
  ggplot(aes(x=paste0(clusters, "(", DeltaMZ , ")"),
             y=max_chemdir, fill=condition)) +
  xlab("Neighboring nodes")+ ylab("ChemProp score")+
  theme_pubr()+
  scale_fill_manual(values = c("#39B54A", "#283891", "#EF3E36"))+
  scale_color_manual(values = c("#39B54A", "#283891", "#EF3E36"))+
  geom_bar(stat='identity', alpha=0.8)+
  facet_grid(ComponentIndex ~ ., scales = "free_y", space="free", switch = "y") + 
  theme(axis.text.y = element_text(size=7),
        panel.spacing.y = unit(0, "lines")) +
  scale_x_discrete(position = "top")+
  coord_flip()
ggsave(paste0(fig_pre,  "main_chemprop_barplot_HEgrade_hep.pdf"), width = 3.5, height = 4)
ggsave(paste0(fig_pre,  "main_chemprop_barplot_HEgrade_hep.png"), width = 3.5, height = 4)

prepost_for_pair <- function(dat, usemedian = T){
  
  if(usemedian == TRUE) {
    ggplot(dat, aes(x=Sample, y=peak_area, 
                    color=as.factor(as.character(row.ID)), group=row.ID)) +
      theme_bw()+
      # scale_fill_manual(values = c("#39B54A", "#283891", "#EF3E36"))+
      #  scale_color_manual(values = c("#39B54A", "#283891", "#EF3E36"))+
      stat_summary(fun = median) +
      stat_summary(fun = median, geom = "line")+
      guides(color=guide_legend(title="Metabolite"))-> p
    return(p)
  }
  else {
    
    ggplot(dat, aes(x=Sample, y=peak_area, 
                    color=as.factor(as.character(row.ID)), group=row.ID)) +
      theme_bw()+
      #scale_fill_manual(values = c("#39B54A", "#283891", "#EF3E36"))+
      #scale_color_manual(values = c("#39B54A", "#283891", "#EF3E36"))+
      stat_summary(fun = mean) +
      stat_summary(fun = mean, geom = "line")+
      guides(color=guide_legend(title="Metabolite"))-> p
    return(p)
  }
}


# Figure 5d - Pre to Post median 

# gpc
subset(lqtab, blood_origin == "hepatic" & row.ID %in% c(33, 190)) %>%
  prepost_for_pair() +
  theme_pubr()+
  scale_color_manual(values = c("#E69F00",  "#999999")) +
  theme(legend.position = "right")
ggsave(paste0(fig_pre,  "median_hep_33_190.pdf"), width = 4, height = 2)

# bilirubins
subset(lqtab, blood_origin == "hepatic" & row.ID %in% c(733, 219)) %>%
  prepost_for_pair() +
  theme_pubr()+
  scale_color_manual(values = c("#E69F00",  "#999999")) +
  theme(legend.position = "right")
ggsave(paste0(fig_pre,  "median_hep_733_219.pdf"), width = 4, height = 2)


# Figure 5f/g - metabolite lineplot by grade
subset(lqtab, blood_origin == "hepatic" & row.ID %in% c(527, 733))%>%
  prepost_for_pair() +
  theme_pubr()+
  scale_color_manual(values = c("#999999", "#E69F00")) +
  facet_grid(~Worst_PostTIPS_HE_mod)+
  theme(legend.position = "bottom")
ggsave(paste0(fig_pre,  "median_byHEgrade_hep_527_733.pdf"), width = 5, height = 3)

subset(lqtab, blood_origin == "hepatic" & row.ID %in% c(29, 406)) %>%
  prepost_for_pair() +
  theme_pubr()+
  scale_color_manual(values = c("#E69F00", "#999999")) +
  facet_grid(~Worst_PostTIPS_HE_mod)+
  theme(legend.position = "bottom")
ggsave(paste0(fig_pre,  "median_byHEgrade_hep_29_406.pdf"), width = 5, height = 3)

# Scores for each pair by grading
subset(grad, group %in% c("hep") & clusters %in% "29_406") 
subset(grad, group %in% c("hep") & clusters %in% "527_733") 

# ---------------------------------------------#
#           Response to reviewers
# ---------------------------------------------#

source("./code/rtr_figures.R")
