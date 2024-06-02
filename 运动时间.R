library(TwoSampleMR)
library(readxl)
# online
setwd('E:\\运动时间\\暴露')
files <- list.files('E:\\运动时间\\暴露')[1:13]
disease_id <- read.csv("C:\\Users\\86182\\.spyder-py3\\mr_ieu_id\\disease_id.csv",header = F)[[1]]
disease_ref <- read.csv("C:\\Users\\86182\\.spyder-py3\\mr_ieu_id\\disease_id.csv",header = F)
for (d in disease_id[192:518]) {
  if (strsplit(d,'[-]')[[1]][1] == 'ebi'){
    next
  }
  for (f in files) {
    exposure_online <- read_excel(f)
    outcome_online <- extract_outcome_data(snps = exposure_online$SNP,outcomes = d)
    dat<-try(harmonise_data(exposure_dat = exposure_online,outcome_dat = outcome_online))
    if ("try-error" %in% class(dat)){
      next
    }
    dat$exposure <- strsplit(f,'[_]')[[1]][1]
    dat$outcome <- disease_ref[d==disease_ref$V1,3]
    mr <- mr(dat, method_list = c('mr_ivw'))
    write.table(mr,file = 'result.csv',append = T,sep = ',',row.names = F,col.names = F)
  }
  Sys.sleep(30)
}




for (diease in files[1:2]) {
  c <- fread(diease)
  for (exposure in exposure_name) {
    k <- data.frame(id=1:10)
    exposure_file <- read_excel(paste0('E:\\运动时间\\暴露\\',exposure))
    d<-merge(exposure_file,c,by = "SNP")
    write.csv(d,file = "C:\\Users\\86182\\AppData\\Local\\R\\win-library\\4.3\\TwoSampleMR\\outcome.csv")
    out<-system.file("outcome.csv",package = "TwoSampleMR")
    outcome_dat<-read_outcome_data(snps = exposure_file$SNP,
                                   filename = out,
                                   sep = ",",
                                   snp_col = "SNP",
                                   beta_col = "beta",
                                   se_col = "se",
                                   eaf_col = "eaf",
                                   effect_allele_col = "effect allele",
                                   other_allele_col = "other allele",
                                   pval_col = "pval")
    dat<-harmonise_data(exposure_dat = exposure_file,outcome_dat = outcome_dat)
    dat$exposure <- strsplit(exposure,'[_]')[[1]][1]
    dat$outcome <- strsplit(diease,'[.]')[[1]][1]
    k$exposure <- dat$exposure[1]
    k$outcome <- dat$outcome[1]
    mr <- mr(dat)
    k$p <- mr[mr$method=='Inverse variance weighted',9]
    k$nsnp <- mr[mr$method=='Inverse variance weighted',6]
    mr_odd <- generate_odds_ratios(mr_res = mr)
    k$OR <- mr_odd$or[mr_odd$method=='Inverse variance weighted']
    k$OR_L <- mr_odd$or_lci95[mr_odd$method=='Inverse variance weighted']
    k$OR_u <- mr_odd$or_uci95[mr_odd$method=='Inverse variance weighted']
    writexl::write_xlsx(dat,paste0('E:\\运动时间\\结果3\\','dat_',strsplit(exposure,'[_]')[[1]][1],'_',strsplit(diease,'[.]')[[1]][1],'.xlsx'))
    writexl::write_xlsx(mr_odd,paste0('E:\\运动时间\\结果3\\','mr_',strsplit(exposure,'[_]')[[1]][1],'_',strsplit(diease,'[.]')[[1]][1],'.xlsx'))
    presso <- mr_presso(BetaExposure = 'beta.exposure',BetaOutcome = 'beta.outcome',SdExposure = 'se.exposure',SdOutcome = 'se.outcome',dat,OUTLIERtest = T, DISTORTIONtest = T)
    k$presso <- presso[2]$`MR-PRESSO results`$`Global Test`$Pvalue
    mr_scatter_plot(mr,dat)
    ggsave(paste0(paste0('E:\\运动时间\\结果3\\','scatter_',strsplit(exposure,'[_]')[[1]][1],'_',strsplit(diease,'[.]')[[1]][1],'.pdf')),width = 5,height = 4)
    mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
    ggsave(paste0(paste0('E:\\运动时间\\结果3\\','funnel_',strsplit(exposure,'[_]')[[1]][1],'_',strsplit(diease,'[.]')[[1]][1],'.pdf')),width = 5,height = 4)
    k$pleiotropy <- mr_pleiotropy_test(dat)$pval
    k$heterogeneity <- mr_heterogeneity(dat)$Q_pval[2]
    mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
    ggsave(paste0(paste0('E:\\运动时间\\结果3\\','leaveoneout_',strsplit(exposure,'[_]')[[1]][1],'_',strsplit(diease,'[.]')[[1]][1],'.pdf')),width = 5,height = 4)
    l <- rbind(l,k[1,])
  }
}
writexl::write_xlsx(l,'E:\\运动时间\\mr_result.xlsx')




