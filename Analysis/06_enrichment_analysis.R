#06_enrichment_analysis.R

rm(list=ls())

library(dplyr)
library(ggplot2)
#1.SNP-based PheWAS----
#41 lead SNP
SNP_lead=read.table('./GWAS/pheWAS/26_Topic_lead_SNP.txt')
SNP=list()
for(i in SNP_lead$V1){
  if(file.exists(paste0('./GWAS/pheWAS/',i,'.csv'))){
    SNP[[i]]=read.csv(paste0('./GWAS/pheWAS/',i,'.csv'),header = T)
    SNP[[i]]$SNP=i
  }else{
    message(i, 'no file!')
  }
}

all=read.csv('./gwasATLAS/4756_GWAS_202309.csv')
unique(all$Domain)
unique(all$Trait)
num_taxa=all%>%group_by(Domain)%>%summarise(n_domain=length(atlas.ID),BgRatio=n_domain/4756)%>%ungroup()

#fisher exact test
results=list()
for(i in names(SNP)){
  print(paste0('=======',i,'======='))
  our=nrow(SNP[[i]])
  
  results[[i]]=data.frame(SNP=i,Domain=num_taxa$Domain,GeneRatio=NA,BgRatio=num_taxa$BgRatio,
                          enrichment_fold=NA,pvalue=NA,k=NA,n=NA,M=num_taxa$n_domain,N=4756)
  
  for(j in 1:nrow(num_taxa)){
    #each domain
    background=num_taxa$n_domain[j]
    temp2=subset(SNP[[i]],Domain==num_taxa$Domain[j])#交集
    
    GeneRatio=nrow(temp2)/our
    enrichment_fold=round(GeneRatio/num_taxa$BgRatio[j],2)
    
    p=phyper(nrow(temp2)-1, background, 4756-background, our, lower.tail = F)
    message(num_taxa$Domain[j],' P: ',round(p,digits = 4))
    results[[i]][j,c('GeneRatio','enrichment_fold','pvalue','k','n')]=c(GeneRatio,enrichment_fold,round(p,digits = 4),nrow(temp2),our)
  }
}

#FDR-corrected P
kk=do.call(rbind,results)
kk$padj=p.adjust(kk$pvalue,method ='fdr')
kk$Sig=case_when(kk$padj<0.001~'***',
                 kk$padj<0.01~'**',
                 kk$padj<0.05~'*',
                 kk$padj<0.1~'.',
                 TRUE~'')


kk_ef=reshape2::dcast(kk,SNP~Domain,value.var = 'enrichment_fold')
kk_p=reshape2::dcast(kk,SNP~Domain,value.var = 'padj')
colnames(kk_p)[2:30]

colSums(kk_p[,2:30]<0.05)%>%sort()
rownames(kk_p)=kk_p$SNP
rowSums(kk_p[,2:30]<0.05)%>%sort()

remain_id=names(which(colSums(kk_p[,2:30]<0.05)>0))#12
remain_SNP_id=kk_p$SNP[which(rowSums(kk_p[,2:30]<0.05)>0)]#22
kk=subset(kk,Domain%in%remain_id&SNP%in%remain_SNP_id)

#visualization
summary(kk$enrichment_fold)

rownames(kk_ef)=kk_ef$SNP
kk_ef=kk_ef[remain_SNP_id,remain_id]

celltype_d <- dist(kk_ef)
celltype_fit <- hclust(celltype_d, method="ward.D")
order_celltypes = rownames(kk_ef)[celltype_fit$order]

genera_d <- dist(t(kk_ef))
genera_fit <- hclust(genera_d, method="ward.D")
order_genera = colnames(kk_ef)[genera_fit$order]

P1=ggplot(kk,aes(x=Domain,y=SNP)) +
  geom_tile(aes(fill=enrichment_fold),color='white') +
  geom_text(aes(label=Sig),size=5) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size=14,color='black',angle=90,hjust=1),
    axis.text.y = element_text(size=14,color='black'),
    axis.title = element_text(size=16,color='black'),
    strip.text = element_text(size=14),
    plot.title = element_text(hjust=0.5,size=18),
    axis.ticks = element_blank()
  ) + 
  scale_y_discrete(limits=order_celltypes) +
  scale_x_discrete(limits=order_genera)+
  labs(y="",title='') + scale_fill_gradientn(colours=c('white','pink','red')) + coord_fixed()

#2.gene-based enrichment_analysis----
group_col=c("#982b2b","#db6968","#EDB3B3","#D0E5D0","#459943")
names(group_col)=c('AD','MCI','SCD','SCS','NC')
gene=read.table('./02_figure/11_FUMA_All/174PC_gene_LD/geneIDs.txt',header = T)#174

kk0=read.csv(file='./02_figure/11_FUMA_All/174PC_gene_LD/GS.txt',header = T,sep = '\t')
table(kk0$Category)

our=174
kk0$Jaccard=kk0$N_overlap/(kk0$N_genes+our-kk0$N_overlap)
summary(kk0$Jaccard)#max 0.1

#biological functions----
kk=subset(kk0,Category%in%c('GO_bp','Reactome'))#20
kk$GeneSet=gsub('GOBP_|REACTOME_','',kk$GeneSet)

kk$GeneSet=tolower(gsub("_", " ",kk$GeneSet))
kk$GeneSet[kk$Category=='Reactome'&kk$GeneSet=="sensory perception"]="sensory perception2"
kk$GeneSet=factor(kk$GeneSet,levels = kk$GeneSet)#"sensory perception" 有两个重复
max(kk$Jaccard)

P1=ggplot(kk, aes(y = -log10(adjP), x = factor(GeneSet,levels = rev(GeneSet)),fill=Jaccard)) + 
  geom_bar(stat = "identity")+
  facet_grid(Category~., scales="free", space="free")+coord_flip()+
  theme_classic() + scale_fill_gradientn(colours=c(c("#3288bd", "#d53e4f")),limits = c(0,0.07)) +
  labs(x = "", y = "-log10(P)")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")

#GWAS_catelog----
kk=subset(kk0,Category%in%c('GWAScatalog'))#52
max(kk$Jaccard)
kk=kk[1:21,]
kk$GeneSet=factor(kk$GeneSet,levels = kk$GeneSet)

P2=ggplot(kk, aes(y = -log10(adjP), x = factor(GeneSet,levels = rev(GeneSet)),fill=Jaccard)) + 
  geom_bar(stat = "identity")+
  coord_flip()+
  theme_classic() + scale_fill_gradientn(colours=c(c("#3288bd", "#d53e4f")),limits = c(0,0.07)) +
  labs(x = "", y = "-log10(P)")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")

P=ggpubr::ggarrange(P1,P2, ncol=2, nrow=1,common.legend = T,legend = 'right') 
