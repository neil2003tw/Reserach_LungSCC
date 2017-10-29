library(gplots)
setwd('/Users/NeilWu/Github/Research_LungSCC_prognosis/')

z_trans<-function(l){
  x<-l[!is.na(l)]
  y<-(x-mean(x))/sd(x)
  l[!is.na(l)]<-y
  return(l)
}

table_cox<-read.table('Cox_list_AD_total.txt',sep='\t')
table_coxb<-read.table('Coxb_list_AD_total.txt',sep='\t')
head(table_cox)
colSums(is.na(table_cox[,-1]))
28833-colSums(is.na(table_cox[,-1])) #SC
#31682-colSums(is.na(table_cox[,-1])) #AD

table_coxsig<-table_cox
for(i in 2:dim(table_coxsig)[2]){
  table_coxsig[is.na(table_coxsig[,i]),i]<-1
  table_coxsig[,i]<-as.numeric(table_coxsig[,i]<0.05)
  #table_coxsig[,i]<-as.numeric(table_coxsig[,i]<0.02)
}
table_coxsig$TOTAL<-apply(table_coxsig[,-1],1,function(x) paste(x,collapse = ''))
colSums(table_coxsig[,-c(1,11)])
colSums(table_coxsig[,-c(1,9)]) #AD
colSums(table_coxsig[,-c(1,11)])/(28833-colSums(is.na(table_cox[,-1]))) #SC
colSums(table_coxsig[,-c(1,9)])/(31682-colSums(is.na(table_cox[,-1]))) #AD

#table_coxsigsub<-table_coxsig[,c(3:5,7,9)]
#table_coxsigsub<-table_cox[,-c(1,11)]
table_coxsigsub<-table_coxsig[,-1][,c(3,4,5,7)]#[,idx1] ##SC
#table_coxsigsub<-table_coxsig[,-1][,c(4:7)]
table_coxbsub<-table_coxb[,-1][,c(3,4,5,7)]#[,idx1] ##SC
#table_coxbsub<-table_coxb[,-1][,c(4:7)]
rownames(table_coxsigsub)<-table_cox[,1]
rownames(table_coxbsub)<-table_coxb[,1]
table_coxbsub<-table_coxbsub[which(rowSums(table_coxsigsub)>2),] #SC
#table_coxbsub<-table_coxbsub[rowSums(table_coxsigsub)>2,] #AD
idx_wrong<-match('day_to_event_death',rownames(table_coxbsub))
table_coxbsub<-table_coxbsub[-idx_wrong,]
#####table_coxsigsubb<-log2(table_coxsigsub)
table_coxbsubb<- table_coxbsub
table_coxbsubb<-data.frame(apply(table_coxbsub,2,z_trans))
#colnames(table_coxbsubb)<-sapply(strsplit(colnames(table_coxbsubb),'cox_beta_'),function(x) x[2]) #AD
colnames(table_coxbsubb)<-sapply(strsplit(colnames(table_coxbsubb),'cox_p_'),function(x) x[2])

max(table_coxbsubb,na.rm = T)
min(table_coxbsubb,na.rm = T)

idx1<-which(0.15>((colSums(table_coxsig[,-c(1,11)]))/(28833-colSums(is.na(table_cox[,-1])))) & 0.05<((colSums(table_coxsig[,-c(1,11)]))/(28833-colSums(is.na(table_cox[,-1])))))
#idx1<-which(0.15>((colSums(table_coxsig[,-c(1,9)]))/(31682-colSums(is.na(table_cox[,-1])))) & 0.05<((colSums(table_coxsig[,-c(1,9)]))/(31682-colSums(is.na(table_cox[,-1])))))

heatmap.2(as.matrix(table_coxbsubb),col = redgreen(100), breaks=seq(-2,2,4/100),
          symm=F,symkey=F,symbreaks=T,cexRow=0.7,cexCol=1,trace="none",srtCol=45)
par(margin(4,3,2,1))

