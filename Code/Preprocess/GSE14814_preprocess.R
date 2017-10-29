setwd('/Users/NeilWu/Github/Research_LungSCC_prognosis/GSE14814_main/')
GSE_data<-read.table('GSE14814_series_matrix.txt',skip = 39,fill = T,sep='\t',stringsAsFactors = F,
                     header = T)

GSE_clinical<-GSE_data[c(15,16,14),]
GSE_clinical[2,]<-grepl('OS status: Dead',GSE_clinical[2,])
GSE_clinical[1,]<-sapply(GSE_clinical[1,],USE.NAMES = F,
                         function(x) strsplit(x,'OS time: ')[[1]][2])
GSE_clinical[1,]<-as.numeric(GSE_clinical[1,])*30
GSE_clinical[3,]<-grepl('Histology type: SQCC',GSE_clinical[3,])

GSE_clinical<-data.frame(t(GSE_clinical[,-1]))
colnames(GSE_clinical)<-c('day_to_event_death','event_death','Subtype_sq')
GSE_clinical[,1]<-as.numeric(as.character(GSE_clinical[,1]))*30

GSE_data<-GSE_data[-c(1:45),]
rownames(GSE_data)<-GSE_data[,1]
GSE_data<-GSE_data[,-1]

GSE_annot<-read.csv('../GSE4573_main/HG-U133A.na35.annot.csv',skip=25,stringsAsFactors = F)
GSE_annot<-GSE_annot[,c('Probe.Set.ID','Gene.Symbol')]
GSE_annot[,2]<-gsub(' ',"",GSE_annot[match(rownames(GSE_data)[-c(1,22285)],GSE_annot[,1]),2])
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
GSE_data_all<-cbind(GSE_data_annotated,
                    GSE_clinical[match(row.names(GSE_data_annotated),row.names(GSE_clinical)),])
GSE_data_merged_sq<-GSE_data_all[ as.character(GSE_data_all$Subtype_sq)=='TRUE',]
GSE_data_merged_sq<-GSE_data_merged_sq[,-dim(GSE_data_merged_sq)[2]]
write.table(GSE_data_merged_sq,'Data_merged_14814.txt',sep='\t',quote = F)

