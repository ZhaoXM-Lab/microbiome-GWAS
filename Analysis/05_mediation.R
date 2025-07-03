#05_mediation.R

rm(list=ls())

library(mediation)
library(dplyr)

meta=read.csv('./01_taxa/all/state_table/00_meta_all.csv')
Topic=read.table('./01_taxa/all/Topic/pheno.txt',header = T)
datas=merge(meta,Topic,by=c('FID','IID'),sort=F)
datas=datas%>%mutate(
  C1.sd = C1/sd(C1, na.rm = T),
  C2.sd = C2/sd(C2, na.rm = T),
  C3.sd = C3/sd(C3, na.rm = T),
  C4.sd = C4/sd(C4, na.rm = T),
  C5.sd = C5/sd(C5, na.rm = T),
  age.sd = Age/sd(Age, na.rm = T),
  gender=relevel(factor(gender),ref='female')
)

GE <- data.table::fread("./01_taxa/all/Topic/test_Topic_207.raw", header=TRUE)
datas=merge(datas,GE[,-c(3:6)],by=c('FID','IID'),sort=F)
rownames(datas)=datas$FID


#2.MMSE/MoCA_B/AD----
results=data.frame(X=NA,M=NA,Y=NA,
                   PathA_beta=NA,PathA_p=NA,
                   PathB_beta=NA,PathB_p=NA,
                   PathC2_beta=NA,PathC2_p=NA,
                   PathC_beta=NA,PathC_p=NA,
                   PathAB_beta=NA,PathAB_p=NA,
                   classify=NA,percentage=NA)

#lead_SNP
SNP_lead=read.table('./01_taxa/all/Topic/26_Topic_lead_SNP.txt')
for(select_cut in SNP_lead$V1){
  select_cut=colnames(GE)[grep(paste0(select_cut,'_'),colnames(GE))]
  print(select_cut)
  
  #MMSE---
  pheno.cov=datas[datas$group%in%c('AD','NC'),c("FID",select_cut,'Topic_8',paste0('C',1:5,'.sd'),'age.sd','gender')]
  a <- lm(Topic_8 ~ .,data=pheno.cov[,-1]) #lm(M~X,df)
  PA=summary(a)
  
  pheno.cov=datas[datas$group%in%c('AD','NC'),c("FID",select_cut,'MMSE','Topic_8',paste0('C',1:5,'.sd'),'age.sd','gender')]#
  b <- lm(MMSE ~., data=pheno.cov[,-1]) #lm(Y~X+M)
  PB=summary(b)#b,c'
  
  c <- lm(MMSE ~., data=pheno.cov[,!colnames(pheno.cov)%in%c("FID", "Topic_8")])
  PC=summary(c)#c
  
  set.seed(123)
  result = mediate(a,b,treat=select_cut,mediator = "Topic_8",boot = T)#默认1000次抽样
  PAB=summary(result)
 
  results=rbind(results,data.frame(X=select_cut,M='ES_Ana',Y='MMSE',
                                   PathA_beta=PA$coefficients[2,'Estimate'],PathA_p=PA$coefficients[2,'Pr(>|t|)'],
                                   PathB_beta=PB$coefficients['Topic_8','Estimate'],PathB_p=PB$coefficients['Topic_8','Pr(>|t|)'],
                                   PathC2_beta=PB$coefficients[2,'Estimate'],PathC2_p=PB$coefficients[2,'Pr(>|t|)'],
                                   PathC_beta=PC$coefficients[2,'Estimate'],PathC_p=PC$coefficients[2,'Pr(>|t|)'],
                                   PathAB_beta=PAB$d.avg,PathAB_p=PAB$d.avg.p,
                                   classify=sign(PAB$d.avg*PB$coefficients[2,'Estimate']),
                                   percentage=ifelse(sign(PAB$d.avg*PB$coefficients[2,'Estimate'])>0,
                                                     round(PAB$d.avg/PAB$tau.coef*100,digits = 3), #ab/c
                                                     round(abs(PAB$d.avg/PAB$z.avg*100),digits = 3)) #ab/c'
  ))
  
  
  #MoCA_B----
  pheno.cov=datas[datas$group%in%c('AD','NC'),c("FID",select_cut,'Topic_8',paste0('C',1:5,'.sd'),'age.sd','gender')]
  a <- lm(Topic_8 ~ .,data=pheno.cov[,-1]) #lm(M~X,df)
  PA=summary(a)
  
  pheno.cov=datas[datas$group%in%c('AD','NC'),c("FID",select_cut,'MoCA_B','Topic_8',paste0('C',1:5,'.sd'),'age.sd','gender')]
  b <- lm(MoCA_B ~., data=pheno.cov[,-1]) #lm(Y~X+M)
  PB=summary(b)#b,c'
  
  c <- lm(MoCA_B ~., data=pheno.cov[,!colnames(pheno.cov)%in%c("FID", "Topic_8")])
  PC=summary(c)#c
  
  set.seed(123)
  result = mediate(a,b,treat=select_cut,mediator = "Topic_8",boot = T)#默认1000次抽样
  PAB=summary(result)
  
 
  results=rbind(results,data.frame(X=select_cut,M='ES_Ana',Y='MoCA_B',
                                   PathA_beta=PA$coefficients[2,'Estimate'],PathA_p=PA$coefficients[2,'Pr(>|t|)'],
                                   PathB_beta=PB$coefficients['Topic_8','Estimate'],PathB_p=PB$coefficients['Topic_8','Pr(>|t|)'],
                                   PathC2_beta=PB$coefficients[2,'Estimate'],PathC2_p=PB$coefficients[2,'Pr(>|t|)'],
                                   PathC_beta=PC$coefficients[2,'Estimate'],PathC_p=PC$coefficients[2,'Pr(>|t|)'],
                                   PathAB_beta=PAB$d.avg,PathAB_p=PAB$d.avg.p,
                                   classify=sign(PAB$d.avg*PB$coefficients[2,'Estimate']),
                                   percentage=ifelse(sign(PAB$d.avg*PB$coefficients[2,'Estimate'])>0,
                                                     round(PAB$d.avg/PAB$tau.coef*100,digits = 3), #ab/c
                                                     round(abs(PAB$d.avg/PAB$z.avg*100),digits = 3)) #ab/c'
  ))
  
  #AD----
  pheno.cov=datas[datas$group%in%c('AD','NC'),c("FID",select_cut,'Topic_8',paste0('C',1:5,'.sd'),'age.sd','gender')]
  a <- lm(Topic_8 ~ .,data=pheno.cov[,-1]) 
  PA=summary(a)
  
  pheno.cov=datas[datas$group%in%c('AD','NC'),c("FID",select_cut,'group','Topic_8',paste0('C',1:5,'.sd'),'age.sd','gender')]#
  pheno.cov$group=relevel(factor(pheno.cov$group),ref='NC')
  b <- glm(group ~., data=pheno.cov[,-1] ,family = binomial)
  PB=summary(b)#b,c'
  
  c <- glm(group ~., data=pheno.cov[,!colnames(pheno.cov)%in%c("FID", "Topic_8")],family = binomial)
  PC=summary(c) 
  
  set.seed(123)
  result = mediate(a,b,treat=select_cut,mediator = "Topic_8",boot = T)#默认1000次抽样
  PAB=summary(result)
  
  results=rbind(results,data.frame(X=select_cut,M='ES_Ana',Y='AD',
                                   PathA_beta=PA$coefficients[2,'Estimate'],PathA_p=PA$coefficients[2,'Pr(>|t|)'],
                                   PathB_beta=PB$coefficients['Topic_8','Estimate'],PathB_p=PB$coefficients['Topic_8','Pr(>|z|)'],
                                   PathC2_beta=PB$coefficients[2,'Estimate'],PathC2_p=PB$coefficients[2,'Pr(>|z|)'],
                                   PathC_beta=PC$coefficients[2,'Estimate'],PathC_p=PC$coefficients[2,'Pr(>|z|)'],
                                   PathAB_beta=PAB$d.avg,PathAB_p=PAB$d.avg.p,
                                   classify=sign(PAB$d.avg*PB$coefficients[2,'Estimate']),
                                   percentage=ifelse(sign(PAB$d.avg*PB$coefficients[2,'Estimate'])>0,
                                                     round(PAB$d.avg/PAB$tau.coef*100,digits = 3), #ab/c
                                                     round(abs(PAB$d.avg/PAB$z.avg*100),digits = 3)) #ab/c'
  ))
}

results=results[-1,]
write.csv(results,'./01_taxa/all/Topic/S14_mediation.csv',row.names = F)