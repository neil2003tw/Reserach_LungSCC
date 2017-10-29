library(survival)
dir_seq<-'Expression-Genes/BI__HT_HG-U133A/Level_3/'
dir_metadata<-'METADATA/BI__HT_HG-U133A/'
dir_clinic<-'../TCGA_main/Biotab/Clinical/Biotab/'

data_srdf<-read.table(paste0(dir_metadata,'broad.mit.edu_LUSC.HT_HG-U133A.sdrf.txt'),
                      sep='\t',stringsAsFactors = F,header=T)
data_followup<-read.table(paste0(dir_clinic,'nationwidechildrens.org_clinical_follow_up_v1.0_lusc.txt'),
                          sep='\t',stringsAsFactors = F,header=T)
data_patient<-read.table(paste0(dir_clinic,'nationwidechildrens.org_clinical_patient_lusc.txt'),
                         sep='\t',stringsAsFactors = F,header=T)

list_gene_array<-grep('.level3.data.txt',dir(dir_seq),value = T)
list_barcode<-data_srdf$Comment..TCGA.Barcode.[match(list_gene_array,
                                                     data_srdf$Derived.Array.Data.Matrix.File.1)]
data_merged<-data.frame(gene_id=NA)

for(i in 1:length(list_gene_array)){
  data_temp<-read.table(paste0(dir_seq,list_gene_array[i]),skip = 1,
                        sep='\t',header=T,stringsAsFactors = F)
  colnames(data_temp)<-c('gene_id',list_barcode[i])
  data_merged<-merge(data_merged,data_temp,by='gene_id',all.y = T)
}

list_gene<-data_merged$gene_id
data_merged<-data.frame(t(data_merged[,-1]))
colnames(data_merged)<-list_gene

data_followup_ordered<-data_followup[match(substr(list_barcode,1,12),data_followup$bcr_patient_barcode),]
data_patient_ordered<-data_patient[match(substr(list_barcode,1,12),data_patient$bcr_patient_barcode),]

data_clinic<-data.frame(barcode=substr(list_barcode,1,12))
data_clinic$event_death<-grepl('Dead',data_followup_ordered$vital_status)
data_clinic$day_to_event_death<-as.numeric(data_followup_ordered$death_days_to)
data_clinic$day_to_event_death[is.na(data_clinic$day_to_event_death)]<-data_followup_ordered$last_contact_days_to[is.na(data_clinic$day_to_event_death)]
data_clinic$day_to_event_death<-as.numeric(data_clinic$day_to_event_death)
data_final<-cbind(data_merged,data_clinic[,-1])

write.table(data_final,'Data_merged_TCarray.txt',sep='\t',quote=F)
