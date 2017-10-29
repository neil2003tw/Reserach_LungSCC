val_file<-dir('Validation_TCGAarray_emp50_gs5_0.4/',full.names = T)
val_file<-grep('Val_',val_file,value = T)

val_template<-read.table(val_file[1])
val_template<-data.frame(gene_id=val_template[,1])

for(i in 1:length(val_file)){
  val_temp<-read.table(val_file[i])
  now_vdata<-strsplit(strsplit(val_file[i],split = '.txt')[[1]],split = 'Val_')[[1]][2]
  val_temp<-val_temp[,c(1,2)]
  colnames(val_temp)[2]<-now_vdata
  val_template<-merge(val_template,val_temp, by.x='gene_id', by.y='gene_id')
}

which.max(rowSums(val_template[,-1]<0.05))
write.table(val_template,'Validation_TCGAarray_emp50_gs5_0.4/Result_table_wAD.txt',sep='\t',quote=F,row.names=F)

val_template<-read.table('Validation_TCGAarray_emp50_gs5_0.4/Result_table_wAD.txt',sep='\t',header=T)
val_template_sig<-val_template
val_template_sig[,-1]<-val_template[,-1]<0.05
val_AD<-val_template_sig[,c(1,1+grep('_AD',colnames(val_template_sig[,-1])))]
val_SC<-val_template_sig[,c(1,1+grep('_AD',colnames(val_template_sig[,-1]),invert = T))]

