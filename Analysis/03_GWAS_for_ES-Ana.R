#03_GWAS_for_ES-Ana.R 

rm(list=ls())

#1.mh plot----
results_log <- data.table::fread("./01_taxa/all/Topic/lm.topic.assoc.linear", head=TRUE)
results_log=results_log[,c('SNP','CHR','BP','P')]
a1=subset(results_log,-log10(P)>2)#60872
a1=a1[,c('SNP','CHR','BP','P')]

dd2=read.table('./01_taxa/all/Topic/clumped_results.clumped',header = T)
#41 SNPs

library(CMplot)
CMplot(a1,plot.type = "m",threshold = c(1e-05,5e-08),threshold.col=c('red','black'),col=c("#4197d8",'grey'),
       threshold.lty = c(1,2),threshold.lwd = c(3,3), amplify = T,
       highlight=dd2$SNP, highlight.col="orange", highlight.pch=19,highlight.cex=1,
       highlight.text=dd2$SNP,
       file.name = 'topic_gwas',file = 'pdf')

#2.meta analysis----
library(dplyr)
results_log <- data.table::fread("./02_SNV_common/all/05_topic/lm.Topic_8.assoc.linear", head=TRUE)
dd2=read.table('./01_taxa/all/Topic/clumped_results.clumped',header = T)
replication=subset(results_log,P<0.1&SNP%in%dd2$SNP)

load(file = './01_taxa/all/state_plot/11_Topic_map_gene.RData')
discovery=map_gene
replication=subset(replication,replication$SNP%in%discovery$SNP)
discovery=discovery[match(replication$SNP,discovery$SNP),]

#association in the same direction
remain_id=c()
for(i in 1:nrow(discovery)){
  if(sign(discovery$BETA[i])==sign(replication$BETA[i])){
    remain_id=c(remain_id,i)
  }
}
discovery=discovery[remain_id,]
replication=replication[remain_id,]

se_summary=data.frame(
  link=discovery$SNP,
  discovery_se=discovery$SE,
  replication_se=replication$SE)

HR_summary=data.frame(
  link=discovery$SNP,
  discovery_HR=discovery$BETA,
  replication_HR= replication$BETA)

data=cbind(se_summary,HR_summary[,-1])

met=function(x) {
  library(metafor)
  y=rma(as.numeric(x[4:5]), sei=as.numeric(x[2:3]))
  y=c(y$b,y$beta,y$se,y$zval,y$pval,y$ci.lb,y$ci.ub,y$tau2,y$I2)
  y
}

#Random-effect model using Restricted maximum likelihood estimator
results=data.frame(t(apply(data,1,met)))
colnames(results)=c("b","beta","se","zval","pval","ci.lb","ci.ub","tau2","I2")
results

out=cbind(discovery[,c('CHR','SNP','BP','A1')],data[,-1],results[,-1])
write.csv(out,file='./01_taxa/all/state_table/S9.csv',row.names = F)

#3.allele----
library(ggplot2)
group_col=c("#982b2b","#db6968","#EDB3B3","#D0E5D0","#459943")
names(group_col)=c('AD','MCI','SCD','SCS','NC')

#discovery----
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

#genotype
GE <- data.table::fread("./01_taxa/all/Topic/test_Topic_207.raw", header=TRUE)
datas=merge(datas,GE[,-c(3:6)],by=c('FID','IID'),sort=F)
colnames(datas)
discovery=datas
discovery$cohort='Discovery'

#Replication----
meta=read.csv('./01_taxa/all/state_table/00_meta_all_v2.csv')
Topic=read.table('./02_SNV_common/all/05_topic/AD2023_topic8.txt',header = T)
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

GE <- data.table::fread("./01_taxa/candidate/meta/test_AD2023_v2.raw", head=TRUE)
datas=merge(datas,GE[,-c(3:6)],by=c('FID','IID'),sort=F)
replication=datas
replication$cohort='Replication'
#3.1 rs910312----
p2=ggplot(data=df0, aes(x=factor(rs910312_T, levels=c(0,1,2),labels=c('CC','TC','TT')), y=Topic_8)) +
  facet_grid(~cohort, scales="free", space="free")+
  geom_jitter(aes(color=group,fill=group),#alpha=0.6,
              position=position_jitterdodge(jitter.width = 0.1, 
                                            jitter.height = 0, 
                                            dodge.width = 0.1)) +
  geom_boxplot(alpha=0.2, width=0.1,
               position=position_dodge(width=0.8),
               size=0.75, outlier.colour = NA) +
  geom_violin(alpha=0.2, width=0.9,
              position=position_dodge(width=0.8),
              size=0.75) +
  scale_color_manual(values = group_col) +
  theme_classic() +
  labs(y='ES-Ana probability', x='rs910312') +
  theme(axis.text = element_text(face = 'bold', size = 10),
        text= element_text(face = 'bold', size = 10),
        legend.position = 'none')

#3.2 rs72327387----
p1=ggplot(data=df0[!is.na(df0$rs72327387_G),], aes(x=factor(rs72327387_G, levels=c(0,1,2),labels = c('GCTTT/GCTTT','G/GCTTT','G/G')), y=Topic_8)) +
  facet_grid(~cohort, scales="free", space="free")+
  geom_jitter(aes(color=group,fill=group),#alpha=0.6,
              position=position_jitterdodge(jitter.width = 0.1, 
                                            jitter.height = 0, 
                                            dodge.width = 0.1)) +
  geom_boxplot(alpha=0.2, width=0.1,
               position=position_dodge(width=0.8),
               size=0.75, outlier.colour = NA) +
  geom_violin(alpha=0.2, width=0.9,
              position=position_dodge(width=0.8),
              size=0.75) +
  scale_color_manual(values = group_col) +
  theme_classic() +
  labs(y='ES-Ana probability', x='rs72327387') +
  theme(axis.text = element_text(face = 'bold', size = 10),
        text= element_text(face = 'bold', size = 10),
        legend.position = 'none')