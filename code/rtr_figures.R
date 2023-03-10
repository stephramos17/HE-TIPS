# rtr_figures.R
# Author: Carolina Dantas
# Date: March 09 2023
# Purpose: To perform analysis and generate figures for reponse to reviewers



# ------------------------------------------------------------ #
#                   Functions
# ------------------------------------------------------------ #
# Check results for 0 vs 2 only
check_results_0_vs_2 <- function(dat) {
  
  subset(dat, componentindex != -1) %>% 
    compare_means(formula = log10_dpeak ~ Worst_PostTIPS_HE_mod ,
                  paired=FALSE,
                  method = "wilcox.test", group.by = "componentindex"
                    ) -> clusters 
  return(subset(clusters, p.adj < 0.1))
  
}


# Answer to reviewer - remove grade 1 and see if results are consistent
# For 0 and 2+ only

# ------------------------------------------------------------ #
#                   Density plot 0 vs 2
# ------------------------------------------------------------ #

# Density plot for peripheral 0 vs 2+
pval_per <- ks.test(subset(per_change,  Worst_PostTIPS_HE_mod == "0", select=log10_dpeak)[[1]],
        subset(per_change,  Worst_PostTIPS_HE_mod == "2+", select=log10_dpeak)[[1]])$p.value

subset(per_change, Worst_PostTIPS_HE_mod %in% c("0", "2+")) %>%
  make_dens_plot() +
  scale_color_manual(values=he_colors[c(1,3)]) + 
  scale_fill_manual(values=he_colors[c(1,3)]) +
  ggtitle(paste("peripheral pval=", 
                format(pval_per, scientific = T)))
ggsave(paste0(fig_pre, "dens_per_deltapeak_rtr.pdf")) # p-value = 1.255e-13


# Density plot for hepatic 0 vs 2+
pval_hep <- ks.test(subset(hep_change,  Worst_PostTIPS_HE_mod == "0", select=log10_dpeak)[[1]],
                    subset(hep_change,  Worst_PostTIPS_HE_mod == "2+", select=log10_dpeak)[[1]])$p.value

subset(hep_change, Worst_PostTIPS_HE_mod %in% c("0", "2+")) %>%
  make_dens_plot() +
  scale_color_manual(values=he_colors[c(1,3)]) + 
  scale_fill_manual(values=he_colors[c(1,3)]) +
  ggtitle(paste("peripheral pval=", 
                format(pval_hep, scientific = T)))
ggsave(paste0(fig_pre, "dens_hep_deltapeak_rtr.pdf")) # p-value = 5.71e-08



# ------------------------------------------------------------ #
#                   Significant clusters 0 vs 2+
# ------------------------------------------------------------ #

# Peripheral
subset(per_change, Worst_PostTIPS_HE_mod %in% c("0", "2+")) %>%
  check_results_0_vs_2()
  

# Hepatic
subset(hep_change, Worst_PostTIPS_HE_mod %in% c("0", "2+")) %>%
  check_results_0_vs_2()


# ------------------------------------------------------------ #
#                   Violin plot 0 vs 2
# ------------------------------------------------------------ #

# Peripheral
per_change %>%
  subset(componentindex %in% c(8)) %>%
  make_violinplot_byGrade()

# Original for peripheral
subset(per_change, componentindex %in% c(8) ) %>%
  ggplot( aes(y=log10_dpeak, 
              x=Worst_PostTIPS_HE_mod, 
              fill=Worst_PostTIPS_HE_mod)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme_bw() +
  geom_point( position = position_jitterdodge(), alpha=0.3) +
  geom_violin(alpha=0.7) +
  xlab("HE grade")+
  scale_fill_manual(values=c("#39B54A", "#283891", "#EF3E36"))+
  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("0", "1"), c("0", "2+"), c("1", "2+"))) 

# Revision for 0 vs 2
subset(per_change, componentindex %in% c(8) & Worst_PostTIPS_HE_mod %in% c("0", "2+")) %>%
  ggplot( aes(y=log10_dpeak, 
                  x=Worst_PostTIPS_HE_mod, 
                  fill=Worst_PostTIPS_HE_mod)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme_bw() +
  geom_point( position = position_jitterdodge(), alpha=0.3) +
  geom_violin(alpha=0.7) +
  xlab("HE grade")+
  scale_fill_manual(values=c("#39B54A",  "#EF3E36"))+
  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("0", "2+"))) 
ggsave(paste0(fig_pre, "per_ba_0_vs_2_rtr.pdf"))

subset(hep_change, componentindex %in% c(8) & Worst_PostTIPS_HE_mod %in% c("0", "2+")) %>%
  ggplot( aes(y=log10_dpeak, 
              x=Worst_PostTIPS_HE_mod, 
              fill=Worst_PostTIPS_HE_mod)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme_bw() +
  geom_point( position = position_jitterdodge(), alpha=0.3) +
  geom_violin(alpha=0.7) +
  xlab("HE grade")+
  scale_fill_manual(values=c("#39B54A",  "#EF3E36"))+
  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("0", "2+"))) 
ggsave(paste0(fig_pre, "hep_ba_0_vs_2_rtr.pdf"))



# ------------------------------------------------------------ #
#                   Specific Bile acids -  0 vs 2+
# ------------------------------------------------------------ #


subset(lqtab, row.ID  %in% bas & Sample == "Post-PIV" & Worst_PostTIPS_HE_mod %in% c("0", "2+")) %>%
  plot_bileacids_violin() +
  facet_wrap(~row.ID) +
  stat_compare_means(method = "wilcox.test")
ggsave(paste0(fig_pre, "posthv_ba_0_vs_2_rtr.pdf"), width = 15, height = 6, units = "cm", dpi = 300)

