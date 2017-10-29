library(survival)
library(Hmisc)

data_list<-grep('Data_merged',dir('./',recursive = T),value = T)
#val_data_loc<-data_list[3]
#val_result_loc<-'Validation_TCGAarray_emp50_gs5_0.4/Val_GSE3141_AD.txt'
in_use_model_cox<-'TCGA_array/result_cox_table.txt'
in_use_model_gset<-'TCGA_array/result_cox_gset_table_emp50_gs5.txt'

#### Predictor extract
Result_coxp<-read.table(in_use_model_cox,stringsAsFactors = F)
Result_gset<-read.table(in_use_model_gset,stringsAsFactors = F,header = T)
hsa_biocarta<-read.table('hsa_pid/hsa_biocarta.csv',sep = ';',stringsAsFactors = F)
hsa_kegg<-read.table('hsa_pid/hsa_kegg.csv',sep=';',stringsAsFactors = F)
hsa_reactome<-read.table('hsa_pid/hsa_reactome.csv',sep=';',stringsAsFactors = F)
hsa_lookup<-data.frame(GSET=c(paste0('Bcarta_',hsa_biocarta$V2),paste0('KEGG_',hsa_kegg$V2),paste0('React_',hsa_reactome$V2)),
                       Gene=c(hsa_biocarta$V3,hsa_kegg$V4,hsa_reactome$V3),stringsAsFactors = F)

sig_gset<-Result_gset[Result_gset$emp_p<0.05,1]
sig_gset_glist<-hsa_lookup[hsa_lookup$GSET %in% sig_gset,]
Result_coxp<-Result_coxp[!is.na(Result_coxp$cox_p),]
#Result_coxp$fdr<-p.adjust(Result_coxp$cox_p,"fdr")
#sig_glist<-sapply(Result_coxp$gene_id[Result_coxp$fdr<0.05],
#       function(x) strsplit(x,'\\.')[[1]][1],USE.NAMES = F)
#predict_marker<-sig_gset_glist[sig_gset_glist$Gene %in% sig_glist,]
predict_marker<-sig_gset_glist

pre_selected_list<-read.table('TCGA_array/result_cox_table.txt')
pre_selected_list<-pre_selected_list$gene_id[ pre_selected_list$cox_p<0.05 ]

for(k in 1:length(data_list)){
  val_data_loc<-data_list[k]
  name_id<-strsplit(data_list[k],'\\/')[[1]][1]
  name_ad<-''
  if(grepl('_AD',data_list[k])){ name_ad<-'_AD'}
  val_result_loc<-paste0('Validation_TCGAarray_emp50_gs5_0.4/Val_',name_id,name_ad,'.txt')
  validate_data <- read.table(val_data_loc,sep='\t',stringsAsFactors = F)
  idx_preselected<-match(pre_selected_list,colnames(validate_data))[!is.na(match(pre_selected_list,colnames(validate_data)))]
  validate_data<-validate_data[,c(idx_preselected,(dim(validate_data)[2]-1):dim(validate_data)[2])]

  if(grepl('TCGA_main',val_data_loc)){
  colnames(validate_data)<-sapply(colnames(validate_data),function(x) strsplit(x,'\\.')[[1]][1])
  }
  validate_markersub <- data.frame(Sample=row.names(validate_data))
  predict_marker$GSET <- factor(predict_marker$GSET)

  for(i in 1:length(levels(predict_marker$GSET))){
    index_gset_raw<-match(predict_marker[as.numeric(predict_marker$GSET)==i,2],
                          colnames(validate_data))
    index_gset<-index_gset_raw[!is.na(index_gset_raw)]
  
    if(length(index_gset_raw)==sum(is.na(index_gset_raw))){
      expr_text<-parse(text = paste0('validate_markersub$Gset',i,'<-NA'))
      eval(expr_text)
      next
    }
    data_sub<-data.frame(validate_data[,index_gset])
    for(j in 1:length(index_gset)){
      data_sub[,j]<-as.numeric(cut2(as.numeric(data_sub[,j]),g = 2))-1
    }
    if(length(index_gset)>1){
      expr_text<-parse(text = paste0('validate_markersub$Gset',i,'<-rowSums(data_sub/length(index_gset))'))
    }else{
      expr_text<-parse(text = paste0('validate_markersub$Gset',i,'<-unlist(data_sub,use.names = F)'))
    }
    eval(expr_text)
  }

  Surv_info<-validate_data[,(dim(validate_data)[2]-1):dim(validate_data)[2]]
  Surv_info$day_to_event_death<-as.numeric(as.character(Surv_info$day_to_event_death))
  validate_markersub<-cbind(validate_markersub[,-1],Surv_info)
  cox.p<-c()
  cox.beta<-c()
  cox.se<-c()
  pb<-txtProgressBar(min=1,max=dim(validate_markersub)[2]-2,style=3)
  for( i in 1:(dim(validate_markersub)[2]-2)){
    if(all(is.na(validate_markersub[,i]))){next}
    cox.summary<-summary(coxph(Surv(time = day_to_event_death, event = event_death) ~ validate_markersub[,i],
                               data = validate_markersub))
    cox.p[i]<-cox.summary$coefficient[5]
    cox.beta[i]<-cox.summary$coefficient[1]
    cox.se[i]<-cox.summary$coefficient[3]
    setTxtProgressBar(pb,i)
  }
  result_gset_table<-data.frame(gene_id=colnames(validate_markersub)[1:(dim(validate_markersub)[2]-2)],
                                cox_p=cox.p,cox_beta=cox.beta,cox_se=cox.se)
  write.table(result_gset_table,val_result_loc,sep='\t',quote=F)
}
