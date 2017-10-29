GSE_data<-read.table('/Users/NeilWu/Github/Research_EGFRmut_survival/GSE31210/GSE31210_series_matrix.txt',
                     skip = 25,fill=T,stringsAsFactors = F)
U133plus2<-read.csv('/Users/NeilWu/Github/Research_lncRNA-AML/HG-U133_Plus_2.na34.annot.csv',stringsAsFactors=F)
idx_clinical<-which(grepl('gene alteration status:',GSE_data[17,]))
GSE_data<-GSE_data[,c(1,idx_clinical)]

relapse_status<-c()
for(i in 1:dim(GSE_data)[2]){
  relapse_status_idx<-grep('relapse: ',GSE_data[1:30,i])
  relapse_status[i]<-!grepl('not',GSE_data[relapse_status_idx,i])
}
relapse_day<-c()
for(i in 1:dim(GSE_data)[2]){
  relapse_day_idx<-grep('days before relapse/censor:',GSE_data[1:35,i])
  relapse_day[i]<-unname(as.numeric(sapply(GSE_data[relapse_day_idx,i], function(x) strsplit(x,'days before relapse/censor: ')[[1]][2])))
}

sample_id<-as.character(GSE_data[49,])
clinical_data<-data.frame(sample_id,event_death=relapse_status,day_to_event_death=relapse_day)  
expression_data<-GSE_data[49:(dim(GSE_data)[1]-1),] ## row49: GSM ID

matched_probe_raw<-as.character(U133plus2$Gene.Symbol[match(expression_data[,1],as.character(U133plus2$Probe.Set.ID))])
matched_probe_done<-sapply(matched_probe_raw,function(x) strsplit(x,' /// ')[[1]][1],USE.NAMES = F)
expression_data[,1]<-matched_probe_done

GSE_data_annotated<-data.frame(t(expression_data),stringsAsFactors = F)
GSE_data_annotated<-cbind(GSE_data_annotated[,-1],clinical_data[,-1])
exp_header<-as.vector(t(GSE_data_annotated[1,]))
colnames(GSE_data_annotated)[-c(dim(GSE_data_annotated)[2]-1:0)]<-exp_header[-c(dim(GSE_data_annotated)[2]-1:0)]
GSE_data_annotated<-GSE_data_annotated[-1,]

write.table(GSE_data_simplified,'GSE31210_main/Data_merged_31210.txt',sep='\t',quote = F)
