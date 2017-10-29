library(survival)
library(Hmisc)
setwd('/Users/NeilWu/Github/Research_LungSCC_prognosis/')
tmpfun <- function(x) as.formula(paste("Surv(day_to_event_death,event_death)",x,sep="~"))

data_list<-grep('Data_merged(?=.*_AD)',dir('./',recursive = T),value = T,perl = T)

for(i in 1:length(data_list)){
  current_data<-read.table(data_list[i],sep='\t',stringsAsFactors = F)
  current_data$day_to_event_death<-as.numeric(current_data$day_to_event_death)
  current_dataname<-strsplit(data_list[i],'/')[[1]][1]
  print(current_dataname)
  cox_table<-sapply(colnames(current_data),
                    function(x) summary(coxph(tmpfun(x), data = current_data))$coefficient[c(5,1,3)])
  row.names(cox_table)<-c('cox_p','cox_beta','cox_se')
  write.table(cox_table,paste0('Cox_list_AD/',current_dataname,'_list.txt'),sep='\t',quote=F)
}


