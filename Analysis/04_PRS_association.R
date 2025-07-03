#04_PRS_association.R

rm(list=ls())

#1.data preparations----
PRS=read.table('./01_taxa/all/Topic/v6.all.score',header = T)
datas=merge(PRS,meta,by=c('FID','IID'),sort = F)
datas=datas%>%mutate(
  C1.sd = C1/sd(C1, na.rm = T),
  C2.sd = C2/sd(C2, na.rm = T),
  C3.sd = C3/sd(C3, na.rm = T),
  C4.sd = C4/sd(C4, na.rm = T),
  C5.sd = C5/sd(C5, na.rm = T),
  age.sd = Age/sd(Age, na.rm = T),
  gender=relevel(factor(gender),ref='female')
)
Topic=read.table('./01_taxa/all/Topic/pheno.txt',header = T)
datas=merge(datas,Topic[,c(1,3)],by='FID',sort=F)

#2.AD vs. NC----
my_PRS=function(data=datas,phe=c('Topic_8','MMSE','MoCA_B')){
  temp <- data.frame(Set=NA, Threshold=NA, R2=NA, R2.adj=NA, P=NA, Coefficient=NA, Standard.Error=NA)
  for (j in phe) {
    k_test=intersect(c("X5e.08","X1e.05","X5e.05","X0.0001","X0.0005","X0.001","X0.005","X0.01","X0.05","X0.1","X0.5","X1"),colnames(data))
    for (i in k_test) {
      pheno.cov=data[,c("FID",i,j,paste0('C',1:5,'.sd'),'age.sd','gender')]
      
      colnames(pheno.cov)[2:3]=c('PRS','Pheno')
      pheno.cov$Pheno=resid(lm(Pheno ~ . , data=pheno.cov[,!colnames(pheno.cov)%in%c("FID", "PRS")]))
      pheno.merge <- pheno.cov[,colnames(pheno.cov)%in%c("FID", "Pheno", "PRS")]
      model=summary(lm(Pheno ~ . ,data=pheno.merge[,c("Pheno", "PRS")]))
      
      temp=rbind(temp,c(j,i,model$r.squared,model$adj.r.squared,
                        model$coefficients[2,4],model$coefficients[2,1],model$coefficients[2,2]))
    }
  }
  
  return(temp)
}


datas2=subset(datas,group%in%c('NC','AD'))
results=my_PRS(data=datas2,phe=c('MMSE','MoCA_B','Topic_8'))
results=results[-1,]
results2=results%>%group_by(Set)%>%arrange(desc(R2))
results2=lapply(unique(results$Set),function(x){
  temp=subset(results,Set==x)
  temp=temp[which.max(temp$R2),]
  return(temp)
})
results2=do.call(rbind,results2)
results2$Set[as.numeric(results2$P)<0.05]
#"MMSE" "MoCA_B"  "Topic_8"
results2
xx=results%>%filter(Threshold=='X0.5')

#3.visualization----
#(1) AD vs. NC boxlpot
library(ggsignif)#test = "wilcox.test",
P1=ggplot(data=datas2,aes(x=group,y=X0.5,color=group))+
  geom_jitter(alpha=0.2,
              position=position_jitterdodge(jitter.width = 0.35, 
                                            jitter.height = 0, 
                                            dodge.width = 0.8))+
  geom_boxplot(alpha=0.2,width=0.45,
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.75)+
  scale_color_manual(values = group_col)+
  theme_classic() + labs(x='',y='PRS_AD')+
  theme(legend.position="none") + 
  theme(axis.text = element_text(face = 'bold',size = 10),
        text= element_text(face = 'bold',size = 10),
        axis.title.x = element_blank())+geom_signif(comparisons = list(c("NC", "AD")),
                                                    y_position=c(2.8),
                                                    map_signif_level=TRUE,color='black')

#(2)linear regression----
colnames(datas2)
P1=ggplot(datas2)+geom_point(aes(x=X0.5,y=Topic_8,color=group))+
  scale_color_manual(values = c("#459943","#982b2b"),
                     limits = c('NC','AD')) +
  labs(title ='beta=-0.022; P=0.003',x='PRS_AD',y='ES_Ana')+
  stat_smooth(mapping = aes(x=X0.5,y=Topic_8) ,method=lm)+ theme_classic()+
  theme(text=element_text(face='bold', size=10),
        axis.text = element_text(face = 'bold',size = 10),
        plot.title=element_text(hjust=0.5))

P2=ggplot(datas2)+geom_point(aes(x=X0.5,y=MMSE,color=group))+
  scale_color_manual(values = c("#459943","#982b2b"),
                     limits = c('NC','AD')) +
  labs(title ='beta=-2.187; P=2.768e-05',x='PRS_AD',y='MMSE')+
  stat_smooth(mapping = aes(x=X0.5,y=MMSE) ,method=lm)+ theme_classic()+
  theme(text=element_text(face='bold', size=10),
        axis.text = element_text(face = 'bold',size = 10),
        plot.title=element_text(hjust=0.5))

P=cowplot::plot_grid(P1,P2,nrow = 2,align = "vh")
