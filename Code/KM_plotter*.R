library(data.table)
source('Code/myggsurv.R')

data1<-read.table('GSE14814_main/Data_merged_14814_Chemo.txt',sep='\t',header = T)
data2<-read.table('GSE14814_main/Data_merged_14814_Chemofree.txt',sep='\t',header = T)
data_merged<-rbind(data1,data2)

data1$ARF3_g<-cut2(data1$ARF3,g = 2)

Surv_frame<-data.frame(time=c(data1[,13517],data2[,13517]),
                       event=c(data1[,13518],data2[,13518]),
                       drug=c( as.character(data1$ARF3_g),rep('untreated',nrow(data2))))
                       #drug=c(rep(1,nrow(data1)),rep(0,nrow(data2))))
model_clinic  <-with(TCGA_sub, Surv(time = day_to_event_death/365, event = event_death) ~ Risk) 
#model_clinic <- with(Surv_frame, Surv(time, event) ~ drug)
g_dose <- myggsurv(
  survfit(model_clinic),
  xlab='Time stay in the experiement (years)',
  ylab='Survival probability',
  title='TCGA') + theme_classic(16)

logr<-survdiff(model_clinic)
1 - pchisq(logr$chisq, 1)
