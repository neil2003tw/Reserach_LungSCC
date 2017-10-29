setwd('/Users/NeilWu/Github/Research_LungSCC_prognosis/')
library(survival)
list_dir<-grep('Data_merge',dir(recursive = T),value = T)
list_dir<-grep('Chemo|free|drug',list_dir,invert = T,value = T)
list_dir_AD<-grep('_AD',list_dir,value = T)
list_dir_SC<-grep('_AD',list_dir,value = T,invert = T)

list_dir_U<-list_dir_AD
name_data<-sapply(list_dir_U, function(x) strsplit(x,'\\/')[[1]][1],USE.NAMES = F)
data_des<-data.frame(row.names=name_data)
data_des$SEMA_stat<-NA
data_des$coxp<-NA
data_des$coxb<-NA
data_des$coxse<-NA
for(i in 1:length(list_dir_U)){
  data<-read.table(list_dir_U[i],sep='\t',stringsAsFactors = F)
  idx_SEMA<-grep('SEMA6A',colnames(data))
  if(length(idx_SEMA)==0){
    data_des[name_data[i],'SEMA_stat']<-'None'
    next
  }
  data_des[name_data[i],'SEMA_stat']<-length(idx_SEMA)
  data_sub<-data[,c(idx_SEMA,(ncol(data)-1):ncol(data))]
  result_cox<-coxph(Surv(as.numeric(day_to_event_death),event_death) ~ data_sub[,1], data = data_sub)
  data_des[name_data[i],2:4]<-summary(result_cox)$coefficient[c(5,1,3)]
}

data_des_AD<-data_des

