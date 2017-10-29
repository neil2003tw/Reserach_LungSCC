setwd('/Users/NeilWu/Github/Research_LungSCC_prognosis/TCGA_main/')
dir_seq<-'mRNA_sub/'
dir_metadata<-'METADATA/'
dir_clinic<-'Biotab/Clinical/Biotab/'

data_srdf<-read.table(paste0(dir_metadata,'unc.edu_LUSC.IlluminaHiSeq_RNASeqV2.1.11.0.sdrf.txt'),
                      sep='\t',stringsAsFactors = F,header=T)
data_followup<-read.table(paste0(dir_clinic,'nationwidechildrens.org_clinical_follow_up_v1.0_lusc.txt'),
                          sep='\t',stringsAsFactors = F,header=T)
data_patient<-read.table(paste0(dir_clinic,'nationwidechildrens.org_clinical_patient_lusc.txt'),
                         sep='\t',stringsAsFactors = F,header=T)
data_drug<-read.table(paste0(dir_clinic,'nationwidechildrens.org_clinical_drug_lusc.txt'),
                         sep='\t',stringsAsFactors = F,header=T)
data_drug<-data_drug[-c(1,2),]
data_drug<-aggregate(pharmaceutical_therapy_drug_name ~ bcr_patient_barcode,data = data_drug,function(x) paste(x,collapse = ','))

list_gene_normalized<-grep('.rsem.genes.normalized_results',dir(dir_seq),value = T)

list_barcode<-data_srdf$Comment..TCGA.Barcode.[match(list_gene_normalized,data_srdf$Derived.Data.File)]

data_merged<-data.frame(gene_id=NA)

for(i in 1:length(list_gene_normalized)){
  data_temp<-read.table(paste0(dir_seq,list_gene_normalized[i]),
                        sep='\t',header=T,stringsAsFactors = F)
  colnames(data_temp)[2]<-list_barcode[i]
  ##start from here, process ith data (10/05)  
  data_merged<-merge(data_merged,data_temp,by='gene_id',all.y = T)
}

list_gene<-data_merged$gene_id
data_merged<-data.frame(t(data_merged[,-1]))
colnames(data_merged)<-list_gene

data_followup_ordered<-data_followup[match(substr(list_barcode,1,12),data_followup$bcr_patient_barcode),]
data_patient_ordered<-data_patient[match(substr(list_barcode,1,12),data_patient$bcr_patient_barcode),]
data_drug_ordered<-data_drug[match(substr(list_barcode,1,12),data_drug$bcr_patient_barcode),]

data_clinic<-data.frame(barcode=substr(list_barcode,1,12))
data_clinic$drug<-data_drug_ordered$pharmaceutical_therapy_drug_name
data_clinic$event_death<-grepl('Dead',data_followup_ordered$vital_status)
data_clinic$day_to_event_death<-as.numeric(data_followup_ordered$death_days_to)
data_clinic$day_to_event_death[is.na(data_clinic$day_to_event_death)]<-data_followup_ordered$last_contact_days_to[is.na(data_clinic$day_to_event_death)]

  
data_final<-cbind(data_merged,data_clinic)

data_finall2<-data_final
idx_noise<-c()
for(i in 1:(dim(data_final)[2]-4)){
  if(mean(data_final[,i])<5){
    idx_noise<-c(idx_noise,i)
    next
  }else if(any(data_final[,i]==0)){
    data_finall2[,i]<-data_finall2[,i]+1
  }
  data_finall2[,i]<-log(data_finall2[,i],base = 2)
}
data_finall2<-data_finall2[,-idx_noise]
data_finall2$ACT_status<-grepl('cisplatin|vinorelbin',data_finall2$drug,ignore.case = T)

write.table(data_finall2,'Data_merged_TCGA_drug.txt',sep='\t',quote=F)
