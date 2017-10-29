library(survival)
library(stringr)
setwd('/Users/NeilWu/Github/Research_LungSCC_prognosis/')

list_coxsc<-dir('Cox_list/',full.names = T)
list_coxad<-dir('Cox_list_AD/',full.names = T)

list_cox<-grep('_main_list',list_coxad,value = T)

for(i in 1:length(list_cox)){
  table_coxnow<-read.table(list_cox[i],sep='\t')
  table_coxnow<-data.frame(t(table_coxnow))
  name_now<-strsplit(list_cox[i],'Cox_list_AD\\/\\/|_main_list')[[1]][2]
  colnames(table_coxnow)[2]<-paste0(colnames(table_coxnow)[2],'_',name_now)
  table_coxnow$Gene<-row.names(table_coxnow)
  table_coxnow<-table_coxnow[,c(2,4)]
  row.names(table_coxnow)<-seq(1,dim(table_coxnow)[1])
  if(i==1){
    table_cox<-table_coxnow
  }else{
    table_cox<-merge(table_cox,table_coxnow,by='Gene',all=T)
  }
}

table_coxsig<-table_cox
for(i in 2:dim(table_coxsig)[2]){
  table_coxsig[is.na(table_coxsig[,i]),i]<-1
  table_coxsig[,i]<-as.numeric(table_coxsig[,i]<0.05)
}
table_coxsig$TOTAL<-apply(table_coxsig[,-1],1,function(x) paste(x,collapse = ''))

write.table(table_cox,'Coxb_list_AD_total.txt', sep='\t')

