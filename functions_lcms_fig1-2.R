# Functions ---------------------------------------

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