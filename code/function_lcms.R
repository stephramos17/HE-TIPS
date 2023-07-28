# Functions ---------------------------------------

# 1. General functions --------

load_quant_tab_lcms <- function(inPath = fbmnDir) {
  
  # Input
  #       - inPath is path to FBMN directory with all data from GNPS
  #       - example: inPath = "./data/chemdir_v28_suspect_20220429/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-8068b7cb-view_all_clusters_withID/"
  # Output
  #       - qtab is the quantification table from FBMN
  
  
  fbmn_files <- list.files(inPath, pattern="*tsv$" , recursive=T, full.names = T)
  
  # Read quantification table  
  qtab <- read.csv(list.files(file.path(
    inPath, "quantification_table_reformatted"), full.names = T))  
  return(qtab)
  
}

load_db_tab_lcms <- function(inPath = fbmnDir) {
  
  # Input
  #       - inPath is path to FBMN directory with all data from GNPS
  #       - example: inPath = "./data/chemdir_v28_suspect_20220429/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-8068b7cb-view_all_clusters_withID/"
  # Output
  #       - db is results from FBMN with cluster and network matches
  
  fbmn_files <- list.files(inPath, pattern="*tsv$" , recursive=T, full.names = T)
  
  # read file with cluster/network matches from databases
  dbFile <- grep(pattern="DB_result", x=fbmn_files, value=T)
  db <- read.delim(dbFile, sep = '\t', header = T, na.strings = c("NA", "N/A"))
  return(db)
  
}

load_metadata_lcms <- function(inPath = fbmnDir) {
  
  # Input
  #       - inPath is path to FBMN directory with all data from GNPS
  #       - example: inPath = "./data/chemdir_v28_suspect_20220429/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-8068b7cb-view_all_clusters_withID/"
  # Output
  #       - mdat is metadata available from FBMN related to each sample
  
  # Read metadata
  mdatFile <- file.path(inPath, "metadata_table", "metadata_table-00000.tsv" )
  mdat <- read.delim(mdatFile)
  
  # Edit colnames and levels of factors according to desired plot order
  colnames(mdat) <- gsub("ATTRIBUTE_", "", colnames(mdat))
  mdat$Sample <- factor(mdat$Sample, levels=c("Pre-P", "Pre-HV", "Post-HV", "Post-P", "BD", 
                                              "Pre-HET", "Post-HET", 
                                              "Pre-TIPS Revision","Post-TIPS Revision", "POD-1. TIPS Revision"))
  mdat$blood_procedure <- factor(mdat$blood_procedure,levels=c("pre", "post", "bd", "het") )
  mdat$blood_origin <- factor(mdat$blood_origin,levels=c("peripheral", "hepatic") )
  
  return(mdat)
  
}


get_long_qtab_w_mdat <- function(xqtab=qtab, 
                                 lookup=cdat[c("cluster.index", "componentindex")]) {
  
  # Input
  #       - xqtab is quantificaiton table from FBMN, which can be loaded with load_quant_tab_lcms
  #       - lookup is cluster data information from FBMN
  # Output
  #       - mdat is metadata available from FBMN related to each sample
  
  # Merge tables to add cluster index and annotation
  xqtab <- merge(xqtab, lookup, by.x="row.ID", by.y="cluster.index")
  lqtab <- pivot_longer(xqtab, cols=grep("Peak", colnames(xqtab), value = T),
                        names_to="peak_sample", values_to="peak_area")
  
  mdat$peak_sample <- gsub(" ", ".", mdat$sampleID)
  lqtab <- merge(lqtab, mdat, by="peak_sample", all.x=TRUE )
  lqtab <- lqtab[with(lqtab, order(blood_origin , blood_procedure, -row.ID, pt_id)), ]
  
  return(lqtab)
}


add_log10_lqtab <- function(dat) {
  
  # Input
  #     - dat is quantificaiton table from FBMN
  # Output
  #     - dat_w_log - data with log10(peak_area)
  
  min_peak <- min(subset(dat, peak_area > 0)[,"peak_area"])

  dat %>%
   # mutate(log_var = log10(peak_area + 1)) %>%
    mutate(peak_mod = peak_area + (min_peak*0.5) ) %>%
    mutate(log10 = log10(peak_mod)) ->   dat_w_log 
  
  return(dat_w_log)
}


load_clusterinfo_tab_lcms <- function(inPath=fbmnDir) {
  
  # Input
  #       - inPath is path to FBMN directory with all data from GNPS
  #       - example: inPath = "./data/chemdir_v28_suspect_20220429/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-8068b7cb-view_all_clusters_withID/"
  # Output
  #       - dat is cluster information available from FBMN for features
  
  
  # read file with peak area intensity by file/group for each cluster
  fbmn_files <- list.files(inPath, pattern="*tsv$" , recursive=T, full.names = T)
  inFile <- grep(pattern="clusterinfo_summary", x=fbmn_files, value = T)
  dat <- read.delim(inFile, sep = '\t', header = T, na.strings = c("NA", "N/A")) 
  colnames(dat) <- gsub("ATTRIBUTE_", "", colnames(dat))
  return(dat)
  
}


check_annotation_fbmn <- function(x, dat_db=mdb) {
  
  # Input
  #       x is metabolite id of interest
  #       dat_db is gnps output that has compound information
  # Output
  #       return any matches for x (inlcudes compound name if it exists)

  
  # Check annotation for given enry
  dat <- dat_db[dat_db$X.Scan. %in% x,] # [,c("X.Scan.","Compound_Name" )]
  
  
  return(dat)
}


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n+1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

get_paired_dat_pt <- function(dat) {
  
  # Input
  #       - dat is long quantification table (lqtab) with metadata info
  #             it must contain at least 2 timepoints, like pre-PIV and post-PIV 
  # Output
  #       - pair_dat includes data for all patients that have paired data
  
  # Get all pts w/ pre & post samples
  pt_dat <- unique(dat[,c("pt_id", "blood_origin", "blood_procedure", "Sample")])
  pair_pts <-  unique(pt_dat$pt_id)[unlist(lapply(unique(pt_dat$pt_id), function (x) {
    sub_dat <- subset(pt_dat, pt_id == x)
    all(c("pre", "post") %in% sub_dat$blood_procedure)
  }))]
  
  # Subset data by pt
  pair_dat <- subset(dat, pt_id %in% pair_pts)
  pair_dat <- pair_dat[with(pair_dat, order(blood_procedure, -row.ID, pt_id  )), ]
  
  return(pair_dat)
  
}

add_cluster_subclass_db <- function(dat=db, dat_c=cdat) {
  
  # Output is mdb file, which is db with:
  #           - componentindex from cdat
  #           - merged subclass based on component index
  
  # Add componentindex col from cdat to db
  merged_db <- merge(dat,
                     dat_c[,c("cluster.index", "componentindex")],
                     by.x="X.Scan.", by.y="cluster.index")
  
  # Get subclass2 - from original data table
  merged_db$subclass2 <- merged_db$subclass
  
  # Get subclass3 - from componentindex
  sclass_id <- unique(drop_na(merged_db[,c("subclass", "componentindex")], "subclass"))
  lookup <- subset(sclass_id, componentindex != -1)
  merged_db$subclass3 <- lookup$subclass[match(merged_db$componentindex, lookup$componentindex)]
  
  # Reset subclass column
  merged_db$subclass <- ""
  
  # Populate .......
  
  # NA in both
  all_nas <- apply(is.na(merged_db[,c("subclass2", "subclass3")]), 1, all)
  merged_db[all_nas,]$subclass <- NA
  
  # Same in 2 and 3
  merged_db[which(merged_db[,"subclass2"] == merged_db[,"subclass3"]),]$subclass <- 
    merged_db[which(merged_db[,"subclass2"] == merged_db[,"subclass3"]),]$subclass2
  
  # Named in 2 and NA in 3
  match2 <- !is.na(merged_db$subclass2) & is.na(merged_db$subclass3)
  merged_db[match2,]$subclass <- merged_db[match2,]$subclass2
  
  # Named in 3 and NA in 2
  match3 <- is.na(merged_db$subclass2) & !is.na(merged_db$subclass3)
  merged_db[match3,]$subclass <- merged_db[match3,]$subclass3
  
  merged_db$subclass <- as.factor(merged_db$subclass)
  return(merged_db)
  
}

get_cluster_info_4metabolite <- function(metab, dat) {
  
  # Input
  #       metab -- id of metabolite of interest (example: GUDCA would be metab=515)
  #       dat --  clustersummary data table with info regarding cluster.index, LibraryID, etc
  
  shared_idx <- subset(dat, cluster.index == metab, select = componentindex)[[1]]
  print(shared_idx)
  sub_dat <- subset(dat, componentindex == shared_idx, select = c(cluster.index, LibraryID, componentindex))
  sub_dat <- sub_dat[!sub_dat$componentindex ==  -1,]
  if(nrow(sub_dat) == 0 ){ print("No macthes")}
  return(sub_dat)
  
}


get_metab_subclass  <- function(idx_ids, dat_db) {
  
  # gets subclass of cluster neighbors indexes
  # idx_ids -- cluster.index ids
  print(idx_ids)
  #  metab_subclass <- unique(db[db$X.Scan. %in% idx_ids,]$subclass)
  metab_subclass <- dat_db[dat_db$X.Scan. %in% idx_ids,]
  
  return(metab_subclass)  
  
}

get_id_matches <- function(x, dat_db = mdb, dat_c = cdat) {
  
  # Input if db and cdat tables, along with id to be matched
  # Output is list of matches for each id
  
  if (x %in% dat_db$X.Scan.) {
    mysubclass <- dat_db[dat_db$X.Scan. %in% x,][,"subclass"]
    print(paste0(x, " in DB ", mysubclass))
    print(str(mysubclass))
    return(c("Annotated", x, as.character(mysubclass), "annotated" ))
  }
  
  else {
    id_matches <- get_cluster_info_4metabolite(metab = x, dat = dat_c)
    res <- get_metab_subclass(idx_ids = id_matches$cluster.index, dat_db = dat_db)
    
    if(length(res$subclass) ==1) {
      return(c("Unannotated", x, as.character(unique(res$subclass)), "within annotated cluster" ))
    }
    
    if(length(res$subclass) > 1) {
      res <- drop_na(res, subclass)
      return(c("Unannotated", x, as.character(unique(res$subclass)), "within annotated cluster" ))
    }
    
    else {
      return(c("Unannotated", x, "No matches", "unannotated" ))
    }
  }
}


combine_matched_ids <- function(hits_subclass) {
  
  # Use on output of get_id_matches
  hits_subclass <- do.call(rbind, hits_subclass)
  hits_subclass <- data.frame(hits_subclass)
  colnames(hits_subclass) <- c("annotation", "metab", "subclass", "classification")
  hits_subclass$subclass <- as.factor(hits_subclass$subclass)
  hits_subclass$annotation <- as.factor(hits_subclass$annotation)
  hits_subclass$classification <- as.factor(hits_subclass$classification)
  return(hits_subclass)
}



# Functions for baseline analysis (figure 3) -------------------

get_change_from_baseline <- function(dat) {
  
  # input: long quantification table with filtered paired data
  # output: same table with pre to post change - post data rows only
  
  data_chg <- subset(dat,  blood_procedure %in% c("pre", "post")) %>%
    arrange(row.ID, pt_id, Sample) %>%
    group_by(row.ID, pt_id) %>%
    mutate(change_from_baseline = peak_mod - peak_mod[1L]) %>% # change from baseline
    mutate("d.peak" = peak_mod/peak_mod[1L]) %>% # [1L] -> baseline
    mutate("log10_dpeak" = log10(d.peak)) %>%
    ungroup()
  
  post_change <- subset(data_chg, blood_procedure == "post")
  return(post_change)
}

make_dens_plot <- function(dat){
  
  ggplot(dat, aes(x=log10_dpeak, group=Worst_PostTIPS_HE_mod,
                  fill=Worst_PostTIPS_HE_mod, color=Worst_PostTIPS_HE_mod))+
    theme_pubr()+
    geom_density(adjust=1.5, alpha=.3) +
    scale_color_manual(values=he_colors) + 
    scale_fill_manual(values=he_colors) -> dens_dpeak
  return(dens_dpeak)
}



make_violinplot_byGrade <- function(dat) {
  # input:
  #     - dat: filtered change pre/post data (with log10_dpeak)  
  #             examples input data: per_change; hep_change; all_change
  #             filtering: subset data by componentindex to be plotted
  
  ggplot(dat, aes(y=log10_dpeak, 
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
                       comparisons = list(c("0", "1"), c("0", "2+"), c("1", "2+"))) +
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("Change pre/post")+
    facet_wrap(~componentindex, nrow=1) -> panel_clusters
  
  return(panel_clusters)
  
}

get_sorted_mb_4plot <- function(dat, grade) {
  # input:
  #     - dat is df with baseline change[per_change or hep_change (previous post_change)]
  #     - grade is reference group to order by
  # output is plot
  
  # Order all metabolites by change
  dat  %>%
    group_by(row.ID, Worst_PostTIPS_HE_mod, componentindex) %>%
    summarise(mean_change = mean(log10_dpeak) ) %>%
    arrange(mean_change) %>% 
    as.data.frame() %>%
    mutate(new_order = (1:nrow(.))) -> mc
  
  
  # Create new dataframe with mean rowID by group
  sm_cfb <- dat  %>%
    group_by(row.ID, Worst_PostTIPS_HE_mod, componentindex, blood_origin) %>%
    summarise(mean_change = mean(log10_dpeak))
  
  
  # Get data for group to order by
  s0 <- subset(sm_cfb, Worst_PostTIPS_HE_mod %in% grade )
  
  # Order all by rowID of group "i"
  row_id_0_order <-  s0[order(s0$mean_change),]$row.ID
  sm_cfb$row.ID <- as.factor(as.character(sm_cfb$row.ID))
  sm_cfb$row.ID <- factor(sm_cfb$row.ID, levels=row_id_0_order)
  return(sm_cfb)
}

get_significance_clusters <- function(dat) {
  
  # input: change in baseline dataframe
  # output: significant clusters
  subset(dat, componentindex != -1) %>% 
    compare_means(formula = log10_dpeak ~ Worst_PostTIPS_HE_mod , 
                  method = "kruskal", group.by = "componentindex") -> clusters
  # subset(p < 0.05, select=c(componentindex, p, method, `p.adj`)) -> sig_clusters
  
  #return(sig_clusters)
  return(clusters)
}

run_ks_for_groups <- function(dat) {
  
  c01 <- ks.test(subset(dat,  Worst_PostTIPS_HE_mod == "0", select=log10_dpeak)[[1]],
                 subset(dat,  Worst_PostTIPS_HE_mod == "1", select=log10_dpeak)[[1]])
  c02 <- ks.test(subset(dat,  Worst_PostTIPS_HE_mod == "0", select=log10_dpeak)[[1]],
                 subset(dat,  Worst_PostTIPS_HE_mod == "2+", select=log10_dpeak)[[1]])
  c12 <- ks.test(subset(dat,  Worst_PostTIPS_HE_mod == "1", select=log10_dpeak)[[1]],
                 subset(dat,  Worst_PostTIPS_HE_mod == "2+", select=log10_dpeak)[[1]])
  
  res <- list(`0 vs 1`=c01, `0 vs 2`= c02,`1 vs 2`= c12)
  return(res)
  
}

# Functions for bile acid analysis (figure 4) --------------------

rename_sample_abbreviations <- function(df=lqtab) {
  
  # Rename original sample names from readmission
  # to reconcile similar timepoint between individuals
  
  df$Sample <- recode(df$Sample, `Pre-P` = "Pre-PIV",
                      `Post-P` = "Post-PIV",
                      `Post-TIPS Revision` = "Post-HET",
                      `POD-1. TIPS Revision` = "Post-HET2" )
  
  df <- droplevels(df)
  return(df)
  
}

run_stat_mbt_abundance <- function(dat, test="wilcox") {
  
  res <- compare_means(formula=log10 ~ Worst_PostTIPS_HE_mod, 
                       data = dat, group.by="row.ID", method = test,
                       p.adjust.method = "fdr")
  return(res)
  
} 

plot_bileacids_violin <- function(dat) {
  
  # dat: dataframe w/ log data and HE grade
  dat %>%
    ggplot(aes(y=log10, x=Worst_PostTIPS_HE_mod)) +
    theme_bw() +
    geom_violin() +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    ylab("log10(abundance)")+
    xlab("HE grade")
  
}


# Function for chemdir analysis (figure 5) ---------------

read_chemdir_files <- function(fnames) {
  
  # x is list of files to be read
  purrr::map_df(fnames, function(x) {
    data <- read.delim(x)
    id <- gsub(pattern = ".*CHEMDIR-(.*)-download_data.*", replacement = "\\1", x=x) 
    cbind(file_id = id, data)
  }) -> combined_file
  
  return(combined_file)
  
}

add_new_chemdir_cols <- function(x) {
  
  # x is chemdir data output
  dat <- x
  dat$clusters <- as.factor(paste(dat$CLUSTERID1, dat$CLUSTERID2, sep = "_"))
  dat$abs_dmz <- abs(dat$DeltaMZ)
  dat$chemdir_w_sign <- dat$max_chemdir * dat$sign_chemdir
  return(dat)
}

load_selfloop_table <- function(mydir) {
  
  selfloop_file <- list.files(fbmnDir, pattern="selfloop", 
                              recursive = TRUE, full.names=TRUE)
  print(selfloop_file)
  selfloop_df <- read.delim(selfloop_file, sep="\t", header = TRUE, na.strings = c("", " ", "NA")) 
  
  return(selfloop_df)
} 

read_chemdirjob_mdat <- function(f=file_mdat) {
  
  # Get file metadata info
  fdat <- read.csv(f)[,c("task_id","tool", "group", "timeseries", "HE_status",
                         "condition", "grading_value", "grading_timepoint")]
  fdat$file_id <- as.factor(gsub("^(.{8}).*", "\\1", fdat$task_id))
  print(head(fdat))
  fdat$group <- factor(fdat$group, levels=c("per", "hep"))
  fdat$grading_value <- factor(fdat$grading_value, levels=c("0", "1", "2+"))
  fdat$grading_timepoint <- factor(fdat$grading_timepoint, levels=c("wpreTIPS", "wpostTIPS"))
  print(head(fdat))
  
  return(fdat)
}


plot_geompoint_dmz <- function(dat) {
  
  dat %>%
    ggplot(aes(x=DeltaMZ , y=max_chemdir))+
    theme_pubr() +
    scale_y_continuous(name = "ChemProp Score") + 
    theme(legend.position = "right") +
    ylab("Score") + xlab("Delta m/z") +
    geom_hline(yintercept = 0.5, linetype="dashed", color = "black")+
    geom_point(alpha=0.4, color="grey") +
    facet_wrap(~condition)-> p
  return(p)
}

dmz_edges_colored_plots <- function(dat) {
  # Generate delta m/z plots by subclass
  # or based on cutoff
  
  # Select hits above threshold
  df05 <- subset(dat, max_chemdir > 0.5)
  
  # Find annotations to use
  lapply(unique(df05$ComponentIndex), function(i) {
    
    selection <- subset(mdb, componentindex == i, 
                        select = c(Compound_Name, componentindex, subclass)) 
    if(nrow(selection) >=1) {
      return(i)
    }
    
  }) %>% unlist -> anno_chemdir
  
  
  lapply(unique(df05$ComponentIndex), function(i) {
    
    selection <- subset(mdb, componentindex == i, 
                        select = c(Compound_Name, componentindex, subclass)) 
    if(nrow(selection) >=1) {
      # print(paste0("Found ", i, " in DB" ))
      # print(unique(selection$subclass))
      return(selection$subclass)
    }
    
  })
  
  
  subset(df05, !grepl("adduct", df05$EdgeAnnotation) & 
           abs(max_chemdir) > 1, 
         select= c(EdgeAnnotation, DeltaMZ, abs_dmz, sign_chemdir, max_chemdir )) %>% 
    na.omit %>%
    unique %>%
    arrange(abs_dmz) -> anno_2add
  
  print(anno_2add) 
  
  # Make plots
  dat %>%
    plot_geompoint_dmz() +
    geom_text(
      data = ~ unique(subset(df05, select=c(EdgeAnnotation, DeltaMZ))),
      aes(label = EdgeAnnotation, x = DeltaMZ, y = 1),
      size = 2, angle=90)
  
  # geom_point(data= ~ subset(.x, ComponentIndex %in% anno_chemdir),
  #            aes(color=as.character(ComponentIndex))) + 
  # guides(color=guide_legend(title="Subclass cluster"))
}

make_barplot_edgeannotations <- function(dat ,cutoff = 0.5) {
  
  subset(dat, max_chemdir > cutoff & 
           !grepl("adduct", dat$EdgeAnnotation)) %>%
    drop_na(EdgeAnnotation) %>%
    ggplot(aes(EdgeAnnotation, fill=group))+
    geom_bar(stat="count", alpha=0.7, position="dodge" ) + 
    coord_flip()+
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(paste0("pre vs post >", cutoff, " score")) +
    ylab("Count") +
    scale_fill_manual(values=c("#2E3192", "#EE207C")) + # green: #39B54A"
    facet_wrap(~group)   -> p
  return(p)
}


get_diff_chemdir <- function(dat) {
  # For chemdir > 1, select clusters and
  # get the ones that have the greatest differences when comparing two groups
  # input dat is chemdir df w/ 2 groups
  
  hclusters <- unique(dat[dat$chemdir > 1,]$clusters)
  datf <- dat[dat$clusters %in% hclusters,]
  deltachemdir <- datf[,c("condition", "clusters", "max_chemdir")]
  deltachemdir <- reshape(deltachemdir, idvar ="clusters", timevar="condition", direction="wide")
  deltachemdir$logd <- abs(log(deltachemdir[,2]/deltachemdir[,3]))
  selected <- deltachemdir[deltachemdir$logd > 1,]$clusters
  fdat <- subset(dat, clusters %in% selected & abs_dmz > 1)
  return(fdat)
  
}

chemdir_diff_bygrade_barplot <- function(dat) {
  
  dat %>%
    get_diff_chemdir() %>%
    ggplot(aes(x=paste0(ComponentIndex, "_", clusters, " (DeltaMZ: ", DeltaMZ , ")"),
               y=max_chemdir, fill=condition)) +
    xlab("Neighboring nodes")+
    theme_bw()+
    scale_fill_manual(values = c("#39B54A", "#283891", "#EF3E36"))+
    scale_color_manual(values = c("#39B54A", "#283891", "#EF3E36"))+
    geom_bar(stat='identity', alpha=0.8) +
    coord_flip() -> p
  
  return(p)
  
  
}




# Functions for Figures 1-2 -------------

paired_wilcox<-function(df,compr){
  pr_mtbids.mlms = lapply(df[, unique(rowID)], function (oid) {
    cat('rowID: ', oid, '\n')
    mtb_auc_sub.dt = df[rowID == oid]
    if(compr=="peripheral_vs_hepatic"){
      pval<-wilcox.test(log_counts ~ ATTRIBUTE_blood_origin, data=mtb_auc_sub.dt, paired=TRUE)$p.value
    }
    else{
      pval<-wilcox.test(log_counts ~ ATTRIBUTE_blood_procedure, data=mtb_auc_sub.dt, paired=TRUE)$p.value
    }
  })
  names(pr_mtbids.mlms) = df[, unique(rowID)]
  pval_df=data.table(mtb_id=df[, unique(rowID)],pval=unlist(pr_mtbids.mlms))
  return(pval_df)
}

run_paired_wilcox<-function(df,compr,blood_org,pal_cols){
  if(blood_org=="peripheral"){
    if(compr=="pre_vs_post"){
      mtb_df<-df%>%filter(ATTRIBUTE_blood_origin==blood_org)%>%
        filter(ATTRIBUTE_blood_procedure=="pre"|ATTRIBUTE_blood_procedure=="post")%>%
        filter(!(ATTRIBUTE_pt_id %in% c('1_a','13_a')))%>%
        arrange(ATTRIBUTE_blood_procedure,ATTRIBUTE_pt_id)%>%
        as.data.table()
    }
    else if(compr=="pre_vs_bd"){
      mtb_df<-df%>%filter(ATTRIBUTE_blood_origin==blood_org)%>%
        filter(ATTRIBUTE_blood_procedure=="pre"|ATTRIBUTE_blood_procedure=="bd")%>%
        filter(!(ATTRIBUTE_pt_id %in% c('13_a','19_a')))%>% 
        arrange(ATTRIBUTE_blood_procedure,ATTRIBUTE_pt_id)%>%
        as.data.table()
    }
    else{
      mtb_df<-df%>%filter(ATTRIBUTE_blood_origin==blood_org)%>%
        filter(ATTRIBUTE_blood_procedure=="post"|ATTRIBUTE_blood_procedure=="bd")%>%
        filter(!(ATTRIBUTE_pt_id %in% c('1_a','19_a')))%>% #for post vs bd
        arrange(ATTRIBUTE_blood_procedure,ATTRIBUTE_pt_id)%>%
        as.data.table()
    }
  }
  else if(blood_org=="pre"){
    mtb_df<-df%>%filter(ATTRIBUTE_blood_procedure=="pre")%>%
      filter(!(ATTRIBUTE_pt_id %in% c('1_a','1_b')))%>%
      arrange(ATTRIBUTE_blood_procedure,ATTRIBUTE_pt_id)%>%
      as.data.table()
  }
  else if(blood_org=="post"){
    mtb_df<-df%>%filter(ATTRIBUTE_blood_procedure=="post")%>%
      filter(!(ATTRIBUTE_pt_id %in% c('1_a','1_b','19_a','21_a')))%>%
      arrange(ATTRIBUTE_blood_procedure,ATTRIBUTE_pt_id)%>%
      as.data.table()
  }
  else{
    mtb_df<-df%>%filter(ATTRIBUTE_blood_origin==blood_org)%>%
      filter(ATTRIBUTE_blood_procedure=="pre"|ATTRIBUTE_blood_procedure=="post")%>%
      filter(!(ATTRIBUTE_pt_id %in% c('1_a','13_a','19_a','21_a')))%>%
      arrange(ATTRIBUTE_blood_procedure,ATTRIBUTE_pt_id)%>%
      as.data.table() 
  }
  allpvals_paired<-paired_wilcox(mtb_df,compr)
  padj<-p.adjust(allpvals_paired$pval,method="fdr",n=595)
  allpvals_paired<-allpvals_paired%>%mutate(padj=padj,comparison=compr)
  return(allpvals_paired)
}

paired_plots<-function(df,mtbIDs,compr,blood_org,pval_df,num_cols,pal_cols){
  if(blood_org=="peripheral"){
    if(compr=="pre_vs_post"){
      mtb_df<-df%>%filter(ATTRIBUTE_blood_origin==blood_org)%>%
        filter(ATTRIBUTE_blood_procedure=="pre"|ATTRIBUTE_blood_procedure=="post")%>%
        filter(!(ATTRIBUTE_pt_id %in% c('1_a','13_a')))%>%
        arrange(ATTRIBUTE_blood_procedure,ATTRIBUTE_pt_id)%>%
        as.data.table()
    }
    else if(compr=="pre_vs_bd"){
      mtb_df<-df%>%filter(ATTRIBUTE_blood_origin==blood_org)%>%
        filter(ATTRIBUTE_blood_procedure=="pre"|ATTRIBUTE_blood_procedure=="bd")%>%
        filter(!(ATTRIBUTE_pt_id %in% c('13_a','19_a')))%>% 
        arrange(ATTRIBUTE_blood_procedure,ATTRIBUTE_pt_id)%>%
        as.data.table()
    }
    else{
      mtb_df<-df%>%filter(ATTRIBUTE_blood_origin==blood_org)%>%
        filter(ATTRIBUTE_blood_procedure=="post"|ATTRIBUTE_blood_procedure=="bd")%>%
        filter(!(ATTRIBUTE_pt_id %in% c('1_a','19_a')))%>% #for post vs bd
        arrange(ATTRIBUTE_blood_procedure,ATTRIBUTE_pt_id)%>%
        as.data.table()
    }
  }
  else if(blood_org=="pre"){
    mtb_df<-df%>%filter(ATTRIBUTE_blood_procedure=="pre")%>%
      filter(!(ATTRIBUTE_pt_id %in% c('1_a','1_b')))%>%
      arrange(ATTRIBUTE_blood_procedure,ATTRIBUTE_pt_id)%>%
      as.data.table()
  }
  else if(blood_org=="post"){
    mtb_df<-df%>%filter(ATTRIBUTE_blood_procedure=="post")%>%
      filter(!(ATTRIBUTE_pt_id %in% c('1_a','1_b','19_a','21_a')))%>%
      arrange(ATTRIBUTE_blood_procedure,ATTRIBUTE_pt_id)%>%
      as.data.table()
  }
  else{
    mtb_df<-df%>%filter(ATTRIBUTE_blood_origin==blood_org)%>%
      filter(ATTRIBUTE_blood_procedure=="pre"|ATTRIBUTE_blood_procedure=="post")%>%
      filter(!(ATTRIBUTE_pt_id %in% c('1_a','13_a','19_a','21_a')))%>%
      arrange(ATTRIBUTE_blood_procedure,ATTRIBUTE_pt_id)%>%
      as.data.table() 
  }
  df_sub<-mtb_df%>%filter(rowID %in% mtbIDs)
  df_sub<-df_sub%>%
    mutate(label_name=ifelse(!is.na(Compound_Name),paste("mtb ",rowID," ",Compound_Name,sep=""),
                             paste("mtb ",rowID," m/z",signif(`row m/z`,digits=5))),
           ATTRIBUTE_blood_origin=factor(ATTRIBUTE_blood_origin, levels = c("peripheral", "hepatic")),
           ATTRIBUTE_blood_procedure=factor(ATTRIBUTE_blood_procedure, levels = c("pre", "post","bd")))
  pvals<-pval_df%>%filter(mtb_id %in% mtbIDs)%>%rename(rowID=mtb_id)%>%
    mutate(label=paste("q = ",signif(padj,digits=4)))%>%left_join(.,df_sub,by="rowID")%>%group_by(rowID) %>% 
    filter(row_number()==1)
  p<-ggpaired(df_sub, x ="ATTRIBUTE_blood_procedure", y = "log_counts",
              color="ATTRIBUTE_blood_procedure",line.color = "gray", line.size = 0.4,
              palette = pal_cols)+theme_bw()+theme(legend.position = "none",
                                                   panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(),
                                                   strip.text.x = element_text(size = 6))+
    facet_wrap(~label_name,ncol=num_cols,labeller = label_wrap_gen(multi_line = TRUE),scale="free_y")+
    labs(x="HE visit",y="log10 mtb abun")+geom_text(size= 3, data= pvals,mapping = aes(x = -Inf, y = -Inf, label = label),
                                                    hjust   = -0.1,
                                                    vjust   = -1)
  
  return(p)
}

