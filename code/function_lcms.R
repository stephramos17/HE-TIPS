# Functions ---------------------------------------

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
  
  # read metadata
  mdatFile <- file.path(inPath, "metadata_table", "metadata_table-00000.tsv" )
  mdat <- read.delim(mdatFile)
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
  
  
  dat <- dat_db[dat_db$X.Scan. %in% x,] # [,c("X.Scan.","Compound_Name" )]
  
  
  return(dat)
}


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n+1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

get_paired_dat_pt <- function(dat) {
  
  # dat is long quantification table (lqtab) with metadata info
  
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
  # metab -- id of metabolite of interest
  # dat --  clustersummary data table with info regarding cluster.index, LibraryID, etc
  # ex: gudca
  # cluster index = 515 (metabolite id); component index = 2 (cluster id)
  
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
