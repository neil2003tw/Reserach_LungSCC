library(xlsx)
GSE_data<-read.table('GSE3141_series_matrix.txt',skip = 50,fill = T,sep='\t',stringsAsFactors = F,
                     header = T)
GSE_clinical<-read.table('clinical_data.txt',fill = T, sep='\t',quote = "",header=T,
                         stringsAsFactors = F)
GSE_clinical<-GSE_clinical[,c(1,11,13,16)]
GSE_clinical_sub<-GSE_clinical[ which(GSE_clinical$X98.985=='S'),]

GSE_annot<-read.csv('GPL570-13270.csv',skip=16,stringsAsFactors = F)
GSE_annot<-GSE_annot[,c(1,11)]

rownames(GSE_data)<-GSE_data[,1]
GSE_data<-GSE_data[,-1]

GSE_annot_id<-GSE_annot[!duplicated(GSE_annot[,2]),2]
GSE_annot_id<-GSE_annot_id[!GSE_annot_id=='']

GSE_data_annotated<-GSE_data[1,]
GSE_data_annotated[1,]<-NA
pb<-txtProgressBar(min=1,max=length(GSE_annot_id),style=3)
for(i in 1:length(GSE_annot_id)){
  setTxtProgressBar(pb,i)
  now_gene<-GSE_annot_id[i]
  now_probe_id<-GSE_annot[GSE_annot[,2]==now_gene,1]
  now_sd_list<-apply(GSE_data[now_probe_id,],1,sd)
  now_probe_sd_max<-names(which.max(now_sd_list))
  temp_data<-GSE_data[now_probe_sd_max,]
  rownames(temp_data)<-strsplit(GSE_annot_id[i],'///')[[1]][1]
  GSE_data_annotated<-rbind(GSE_data_annotated,temp_data)
}

GSE_data_annotated<-GSE_data_annotated[-1,]
GSE_data_annotated<-data.frame(t(GSE_data_annotated))
GSE_data_annotated_sub<-GSE_data_annotated[ match(GSE_clinical_sub[,1],
                                                  rownames(GSE_data_annotated)), ]
GSE_clinical_sub[,3]<-gsub('>','',GSE_clinical_sub[,3])
GSE_clinical_sub[,3]<-as.numeric(GSE_clinical_sub[,3])*30
colnames(GSE_clinical_sub)<-c('GSM','Subtype','day_to_event_death','event_death')
GSE_clinical_sub_sub<-GSE_clinical_sub[,c(3,4)]
GSE_data_merged<-cbind(GSE_data_annotated_sub,GSE_clinical_sub_sub)
write.table(GSE_data_merged,'Data_merged_3141.txt',sep='\t',quote = F)



