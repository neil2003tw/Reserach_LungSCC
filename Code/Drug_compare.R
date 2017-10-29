setwd('/Users/NeilWu/Github/Research_LungSCC_prognosis/')
library(caret)
library(survival)

read.coxresult<-function(x){
  y<-read.table(x,header = T,sep='\t',stringsAsFactors = F)
  y<-y[,1:(ncol(y)-2)]
  return(data.frame(t(y)))
  }
z_trans<-function(x){
  y<-(x-mean(x))/sd(x)
  return(y)
 }


GSE14814_ACT<-'Cox_list/GSE14814_main_chlist.txt'
GSE14814_FR<-'Cox_list/GSE14814_main_chfrlist.txt'
TCGA_ACT<-'Cox_list/TCGA_main_chlist.txt'
TCGA_FR<-'Cox_list/TCGA_main_chfrlist.txt'
TCGA_ACTFR<-'Cox_list/TCGA_main_actfrlist.txt'

data_GSE14814_ACT<-read.coxresult(GSE14814_ACT)
data_GSE14814_FR<-read.coxresult(GSE14814_FR)
data_TCGA_ACT<-read.coxresult(TCGA_ACT)
data_TCGA_ACTFR<-read.coxresult(TCGA_ACTFR)
data_TCGA_FR<-read.coxresult(TCGA_FR)

length(which(xor(data_GSE14814_FR$cox_p<0.05,data_GSE14814_ACT$cox_p<0.05)))
sum(data_GSE14814_FR$cox_p<0.05)
sum(data_GSE14814_ACT$cox_p<0.05)

length(which(xor(data_TCGA_ACTFR$cox_p<0.05,data_TCGA_ACT$cox_p<0.05)))
length(which(xor(data_TCGA_FR$cox_p<0.05,data_TCGA_ACT$cox_p<0.05)))
sum(!which(data_TCGA_ACT$cox_p<0.001) %in% which(data_TCGA_FR$cox_p<0.05))
sum(data_TCGA_ACTFR$cox_p<0.05)
sum(data_TCGA_ACT$cox_p<0.05)
sum(data_TCGA_FR$cox_p<0.05)

idx_TCGA<-which(data_TCGA_ACT$cox_p<0.001)[!which(data_TCGA_ACT$cox_p<0.001) %in%
                                           which(data_TCGA_FR$cox_p<0.05)]
!which(data_TCGA_ACT$cox_p<0.001) %in% which(data_TCGA_FR$cox_p<0.05)
TCGA_ACT<-read.table('TCGA_main/Data_merged_TCGA_Chemo.txt',sep='\t',header=T,stringsAsFactors = F)
TCGA_sub<-TCGA_ACT[,c(idx_TCGA,17076,17077)]
colnames(TCGA_sub)<-sapply(strsplit(colnames(TCGA_sub),'\\.'),function(x) x[[1]][1])
TCGA_sub<-TCGA_sub[,-match(c('C22orf28','EVI5L', 'FAM127C','SLC9A9'),colnames(TCGA_sub))]##Additive
TCGA_sub$day_to_event_death<-as.numeric(TCGA_sub$day_to_event_death)
model_cox<-coxph(Surv(time = day_to_event_death, event = event_death) ~ . ,data = TCGA_sub)
TCGA_sub$Risk<-factor(predict(model_cox,TCGA_sub)>0)
levels(TCGA_sub$Risk)<-c('Low-risk','High-risk')
survdiff(Surv(time = day_to_event_death, event = event_death) ~ Risk ,data = TCGA_sub)
TCGA_train<-TCGA_sub[,-c(9,10)]
TCGA_train[,-ncol(TCGA_train)]<-data.frame(apply(TCGA_train[,-ncol(TCGA_train)], 2, z_trans))
inTest<-createDataPartition( TCGA_train$Risk,p=1/4,list=F )
TCGA_train_t<-TCGA_train[-inTest,]
TCGA_train_iv<-TCGA_train[inTest,]
model_svm<-train(Risk ~ .,method='svmLinear',data=TCGA_train_t)
#model_nb<-train(Risk ~ .,method='nb',data=TCGA_train_t)
predictions<-predict(model_svm,TCGA_train_iv)
confusionMatrix(predictions,TCGA_train_iv$Risk)  ## Internal validate

GSE_ACT<-read.table('GSE14814_main/Data_merged_14814_Chemo.txt',sep='\t',header=T)
idx_match<-match(colnames(TCGA_train),colnames(GSE_ACT))
GSE_sub<-GSE_ACT[,c(idx_match[!is.na(idx_match)],13517,13518)]
GSE_sub[,-c(9,10)]<-data.frame(apply(GSE_sub[,-c(9,10)], 2, z_trans))
#GSE_sub$C22orf28<-NA
#GSE_sub$EVI5L<-NA
#GSE_sub$FAM127C<-NA
#GSE_sub$SLC9A9<-NA
predict(model_svm,GSE_sub)
GSE_sub$Risk<-predict(model_svm,GSE_sub)
survdiff(Surv(time = day_to_event_death, event = event_death) ~ Risk ,data = GSE_sub)

GSE14814<-cbind(data_GSE14814_ACT,data_GSE14814_FR)
colnames(GSE14814)[4:6]<-paste0(colnames(GSE14814)[4:6],'_FR')
GSE14814_sel<-GSE14814[GSE14814$cox_p<0.05&GSE14814$cox_p_FR>0.05,]

#GSE14814_selected<-GSE14814[which(data_GSE14814_FR$cox_p<0.05&data_GSE14814_ACT$cox_p<0.05),]
#temp2<-GSE14814_selected[which(GSE14814_selected$cox_beta>0&GSE14814_selected$cox_beta_FR<0),]


TCGA<-cbind(data_TCGA_ACT,data_TCGA_FR)
TCGA<-TCGA[-which(duplicated(sapply(strsplit(row.names(TCGA),'\\.'),function(x) x[[1]][1]))),]
row.names(TCGA)<-sapply(strsplit(row.names(TCGA),'\\.'),function(x) x[[1]][1])
colnames(TCGA)[4:6]<-paste0(colnames(TCGA)[4:6],'_FR')
TCGA_sel<-TCGA[TCGA$cox_p<0.05&TCGA$cox_p_FR>0.05,]



#TCGA_sel<-TCGA[which(TCGA$cox_p<0.05&TCGA$cox_p<0.05),]
#temp4<-TCGA_sel[which(TCGA_sel$cox_beta>0&TCGA_sel$cox_beta_FR<0),]

temp<-TCGA[row.names(TCGA) %in% row.names(GSE14814_sel),]
temp<-temp[temp$cox_p<0.05,]
temp[,4:6]<-GSE14814[match(row.names(temp),row.names(GSE14814)),1:3]
colnames(temp)[4:6]<-c('cox_p_TCGA','cox_beta_TCGA','cox_se_TCGA')

