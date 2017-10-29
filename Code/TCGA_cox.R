Cox_analysis<-function(data_table_loc){
  
  library(survival)
  library(Hmisc)
  library(fmsb)

  data_merged<-read.table(data_table_loc,sep='\t',stringsAsFactors = F)
  data_merged$day_to_event_death<-as.numeric(data_merged$day_to_event_death)
  for(i in 1:(dim(data_merged)[2]-2)){
    data_merged[,i]<-as.numeric(as.character(data_merged[,i]))
  }
  data_merged_sub<-data_merged[!is.na(data_merged$day_to_event_death),]
  #data_merged_sub<-data_merged
  
  ############ Univariate cox for detecting prognosis
  cox.p<-c()
  cox.beta<-c()
  cox.se<-c()
  
  pb<-txtProgressBar(min=1,max=dim(data_merged_sub)[2]-2,style=3)
  for( i in 1:(dim(data_merged_sub)[2]-2)){
    cox.summary<-summary(coxph(Surv(time = day_to_event_death, event = event_death) ~ data_merged_sub[,i],
                               data = data_merged_sub))
    cox.p[i]<-cox.summary$coefficient[5]
    cox.beta[i]<-cox.summary$coefficient[1]
    cox.se[i]<-cox.summary$coefficient[3]
    setTxtProgressBar(pb,i)
  }

  result_table<-data.frame(gene_id=colnames(data_merged_sub)[1:(dim(data_merged_sub)[2]-2)],cox_p=cox.p,cox_beta=cox.beta,cox_se=cox.se)
  write.table(result_table,'result_cox_table.txt',quote=F,sep='\t')
  result_table_sub<-result_table[!is.na(result_table$cox_p),]
  #result_table_sub$fdr<-p.adjust(result_table_sub$cox_p,"fdr")
  ### FDR<0.05 as significant
  result_table_sig<-result_table_sub[ result_table_sub$cox_p<0.1, ]
  
  ########### Midium-seperated classify
  idx_match<-match(result_table_sig$gene_id,colnames(data_merged_sub))
  data_gene_sig<-data_merged_sub[,idx_match[!is.na(idx_match)]]
  data_gene_sig_midium<-data_gene_sig
  for(i in 1:dim(data_gene_sig)[2]){
    data_gene_sig_midium[,i]<-as.numeric(cut2(data_gene_sig[,i],g = 2))-1
  }

  hsa_biocarta<-read.table('../hsa_pid/hsa_biocarta.csv',sep = ';',stringsAsFactors = F)
  hsa_kegg<-read.table('../hsa_pid/hsa_kegg.csv',sep=';',stringsAsFactors = F)
  hsa_reactome<-read.table('../hsa_pid/hsa_reactome.csv',sep=';',stringsAsFactors = F)

  hsa_lookup<-data.frame(GSET=c(paste0('Bcarta_',hsa_biocarta$V2),paste0('KEGG_',hsa_kegg$V2),paste0('React_',hsa_reactome$V2)),
                         Gene=c(hsa_biocarta$V3,hsa_kegg$V4,hsa_reactome$V3))
  result_gene_sig<-sapply(as.character(result_table_sig$gene_id),function(x) strsplit(x,'\\.')[[1]][1],USE.NAMES = F)

  data_gset<-data.frame(Sample=row.names(data_gene_sig))
  hsa_lookup$GSET<-factor(hsa_lookup$GSET)
  Gset_lookup<-levels(hsa_lookup$GSET)

  for(i in 1:length(levels(hsa_lookup$GSET))){
    index_gset_raw<-match(hsa_lookup[as.numeric(hsa_lookup$GSET)==i,2],result_gene_sig)
    #if((sum(is.na(index_gset_raw))/length(index_gset_raw))>0.6){next}
    if(sum(!is.na(index_gset_raw))<2){next}  ##1117 modified, change ratio to minimum count
    else if(length(index_gset_raw)<5){next}
    index_gset<-index_gset_raw[!is.na(index_gset_raw)]
    if(length(index_gset)>1){
      expr_text<-parse(text = paste0('data_gset$Gset',i,'<-rowSums(data_gene_sig_midium[,index_gset])/length(index_gset)'))
    }else{
      expr_text<-parse(text = paste0('data_gset$Gset',i,'<-data_gene_sig_midium[,index_gset]'))
    }
    eval(expr_text)
  }

  ############### Geneset-based Cox regression
  Surv_info<-data_merged_sub[,(dim(data_merged_sub)[2]-1):dim(data_merged_sub)[2]]
  data_gset<-cbind(data_gset[,-1],Surv_info)

  cox.p<-c()
  cox.beta<-c()
  cox.se<-c()

  pb<-txtProgressBar(min=1,max=dim(data_gset)[2]-2,style=3)
  for( i in 1:(dim(data_gset)[2]-3)){
    cox.summary<-summary(coxph(Surv(time = day_to_event_death, event = event_death) ~ data_gset[,i],
                               data = data_gset))
    cox.p[i]<-cox.summary$coefficient[5]
    cox.beta[i]<-cox.summary$coefficient[1]
    cox.se[i]<-cox.summary$coefficient[3]
    setTxtProgressBar(pb,i)
  }

  result_gset_table<-data.frame(gene_id=colnames(data_gset)[1:(dim(data_gset)[2]-3)],cox_p=cox.p,cox_beta=cox.beta,cox_se=cox.se)
  result_gset_table$gene_id<-as.character(result_gset_table$gene_id)
  result_gset_table$gene_id<-Gset_lookup[sapply(result_gset_table$gene_id,function(x) as.numeric(strsplit(x,'Gset')[[1]][2]),USE.NAMES = F)]



######## Geneset based permutation
  result_gset_table$g_count<-NA
  for(i in 1:dim(result_gset_table)[1]){
    index_gset_raw<-match(hsa_lookup[ hsa_lookup$GSET==result_gset_table$gene_id[i],2],result_gene_sig)
    result_gset_table$g_count[i]<-sum(!is.na(index_gset_raw))
    result_gset_table$g_count_total[i]<-table(hsa_lookup$GSET)[ result_gset_table$gene_id[i]]
  }

  result_gset_table$emp_p<-NA
  count_gene<-names(table(result_gset_table$g_count))
  pb<-txtProgressBar(min=1,max=length(count_gene),style=3)
  for(i in 1:length(count_gene)){
    set.seed(12345)
    permute.p<-c()
    setTxtProgressBar(pb,i)
    for(j in 1:50000){  ###Permutation g-set
      idx_permute<-sample(dim(data_gene_sig_midium)[2],as.numeric(count_gene[i]))
      data_temp<-data.frame(data=rowSums(data_gene_sig_midium[,idx_permute]/as.numeric(count_gene[i])))
      data_temp<-cbind(data_temp,Surv_info)
      cox.summary<-summary(coxph(Surv(time = day_to_event_death, event = event_death) ~ data,
                                 data = data_temp))
      permute.p[j]<-cox.summary$coefficient[5]
    }
    idx_current_gcount<-which(result_gset_table$g_count==count_gene[i])
    for(h in idx_current_gcount){
      emp.p<-sum(result_gset_table$cox_p[h]>permute.p)/50000
      result_gset_table$emp_p[h]<-emp.p
    }
  }
  write.table(result_gset_table,'result_cox_gset_table_emp50_gs5.txt',quote=F,sep='\t')
}