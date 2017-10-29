GSE_data<-read.table('GSE4573_series_matrix.txt',skip = 25,fill = T,sep='\t',stringsAsFactors = F)
GSE_clinical<-read.table('GSE4573_clinical_data.txt',fill = T, sep='\t',quote = "",header=T,
                         stringsAsFactors = F)
GSE_annot<-read.csv('HG-U133A.na35.annot.csv',skip=25,stringsAsFactors = F)

GSE_clinical_sub<-GSE_clinical[match(GSE_data[1,],GSE_clinical$RNA_array_SCC_ID),c(6,7)]
GSE_clinical_sub[,2]<-GSE_clinical_sub[,2]*30
colnames(GSE_clinical_sub)<-c('event_death','day_to_event_death')
GSE_clinical_sub<-GSE_clinical_sub[-1,]

GSE_data<-GSE_data[-c(1:33),]
colnames(GSE_data)<-GSE_data[1,]
rownames(GSE_data)<-GSE_data[,1]
GSE_data<-GSE_data[-1,-1]

GSE_annot<-GSE_annot[,c('Probe.Set.ID','Gene.Symbol')]
GSE_annot[,2]<-gsub(' ',"",GSE_annot[match(rownames(GSE_data),GSE_annot[,1]),2])
GSE_annot_id<-GSE_annot[!duplicated(GSE_annot[,2]),2]

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
GSE_data_merged<-cbind(GSE_data_annotated,GSE_clinical_sub)
write.table(GSE_data_merged,'Data_merged_4573.txt',sep='\t',quote = F)



