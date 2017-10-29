gset_sig<-as.character(result_gset_table[result_gset_table$emp_p<0.05,1])
kappa_result<-data.frame(Gene_set1=NA,Gene_set2=NA,Kappa_p=NA,Kappa_estimate=NA,Kappa_judge=NA)
for(i in 1:(length(gset_sig)-1)){
  for(j in (i+1):length(gset_sig)){
    gs1<-as.character(hsa_lookup[hsa_lookup[,1]==gset_sig[i],2])
    gs2<-as.character(hsa_lookup[hsa_lookup[,1]==gset_sig[j],2])
    all_vec<-c(gs1,gs2)
    all_vec<-all_vec[!duplicated(all_vec)]
    gs1<-factor(all_vec %in% gs1,levels = c('TRUE','FALSE'))
    gs2<-factor(all_vec %in% gs2,levels = c('TRUE','FALSE'))
    try(summary.kappa<-Kappa.test(gs1,gs2))
    kappa_temp<-data.frame(Gene_set1=i,Gene_set2=j,
                           Kappa_p=summary.kappa$Result$p.value,
                           Kappa_estimate=summary.kappa$Result$estimate,
                           Kappa_judge=summary.kappa$Judgement)
    kappa_result<-rbind(kappa_result,kappa_temp)
  }
}
