setwd('/Users/NeilWu/Github/Research_LungSCC_prognosis/TCGA_main/')
TCGA_data<-read.table('Data_merged_TCGA_drug.txt',stringsAsFactors = F,header = T,sep='\t')

TCGA_data_dr<-TCGA_data[TCGA_data$ACT_status==1,-c(17076,17077,17080)]
write.table(TCGA_data_dr,'Data_merged_TCGA_Chemo.txt',sep='\t',quote = F)

TCGA_data_actfr<-TCGA_data[TCGA_data$ACT_status==0,-c(17076,17077,17080)]
write.table(TCGA_data_actfr,'Data_merged_TCGA_ACTfree.txt',sep='\t',quote = F)

TCGA_data_fr<-TCGA_data[is.na(TCGA_data$drug),-c(17076,17077,17080)]
write.table(TCGA_data_fr,'Data_merged_TCGA_Chemofree.txt',sep='\t',quote = F)


