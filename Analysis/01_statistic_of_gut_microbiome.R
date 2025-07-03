# 01_statistic_of_gut_microbiome.R

rm(list=ls())

library(dplyr)
library(ggplot2)

group_col=c("#982b2b","#db6968","#EDB3B3","#D0E5D0","#459943")
names(group_col)=c('AD','MCI','SCD','SCS','NC')

#1.meta:group-age-sex----
map=read.csv('./03_table/00_meta_all.csv',header=T)
df_count_sex <- map %>%
  group_by(Age, gender) %>%
  summarise(count = n(), .groups = 'drop')

df_count_group <- map %>%
  group_by(Age, group) %>%
  summarise(count = n(), .groups = 'drop')

p <- ggplot() +
  geom_bar(data = df_count_group, aes(x = Age, y = count, fill = factor(group,levels = c('AD','MCI','SCD','SCS','NC'))), stat = "identity", position = "stack") +
  geom_bar(data = df_count_sex, aes(x = Age, y = -count, fill = gender), stat = "identity", position = "stack") +
  scale_fill_manual(values = c("male" = "#E69F00", "female" = "#CC79A7", group_col)) +
  labs(x = "Age", y = "Number of samples", fill = "Legend") +
  theme_classic()+
  theme(axis.title.y = element_text(vjust = 2), legend.position = "bottom") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")
p

#2.Figure S1----
#(1)APOE----
map=read.csv('./03_table/00_meta_all.csv',header=T)#252
map$APOE_status=factor(map$APOE_status,levels = c('0','1','2'))
map$group=factor(map$group,levels = c('NC','SCS','SCD','MCI','AD'))

df1 = map %>% 
  count(group,APOE_status) %>% 
  group_by(group) %>% 
  mutate(prop=prop.table(n))

P=ggplot(df1, aes(group, prop*100, fill = APOE_status))+
  geom_col() +
  geom_text(aes(label = scales::percent(prop)),
            position = position_stack(vjust = 0.5)) +
  theme_classic()+
  guides(fill=guide_legend(title = expression(paste('APOE ',epsilon,'4 carrier'))))+
  #paste这里就是当paste0用,"APOE ε4 status"
  labs(y='Percentages (%)',x='')+
  theme(axis.text = element_text(face = 'bold',size = 10),
        text= element_text(face = 'bold',size = 10))+scale_fill_brewer(palette="Set2")

#(2)alpha----
map=read.csv('./01_taxa/feat_genus.csv',header=T)
map1=map[,c('group','Richness','Shannon','Simpson')]
map1$group=factor(map1$group,levels = c('NC','SCS','SCD','MCI','AD'))

Significant <- ggpubr::compare_means(Richness ~ group, 
                                     data = map1, p.adjust.method = "fdr",ref.group = 'NC',
                                     symnum.args = list(cutpoints = c(0,0.001,0.01,0.05,0.1,1), 
                                                        symbols = c("***","**","*",'.',"ns")))
P1=ggplot(data=map1,aes(x=group,y=Richness,color=group))+
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
  theme_classic() +
  theme(legend.position="none") + 
  theme(axis.text = element_text(face = 'bold',size = 10),
        text= element_text(face = 'bold',size = 10),
        axis.title.x = element_blank())+
  geom_text(data = Significant,aes(x=group2,y=rep(87,4),label=p.signif),col='black',size=5)

Significant <- ggpubr::compare_means(Shannon ~ group, 
                                     data = map1, p.adjust.method = "fdr",ref.group = 'NC',
                                     symnum.args = list(cutpoints = c(0,0.001,0.01,0.05,0.1,1), 
                                                        symbols = c("***","**","*",'.',"ns")))
P2=ggplot(data=map1,aes(x=group,y=Shannon,color=group))+
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
  theme_classic() +
  theme(legend.position="none") + 
  theme(axis.text = element_text(face = 'bold',size = 10),
        text= element_text(face = 'bold',size = 10),
        axis.title.x = element_blank())+geom_text(data = Significant,aes(x=group2,y=rep(3.5,4),label=p.signif),col='black',size=5)

Significant <- ggpubr::compare_means(Simpson ~ group, 
                                     data = map1, p.adjust.method = "fdr",ref.group = 'NC',
                                     symnum.args = list(cutpoints = c(0,0.001,0.01,0.05,0.1,1), 
                                                        symbols = c("***","**","*",'.',"ns")))
P3=ggplot(data=map1,aes(x=group,y=Simpson,color=group))+
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
  theme_classic() +
  theme(legend.position="none") + 
  theme(axis.text = element_text(face = 'bold',size = 10),
        text= element_text(face = 'bold',size = 10),
        axis.title.x = element_blank())+geom_text(data = Significant,aes(x=group2,y=rep(1,4),label=p.signif),col='black',size=5)

P=cowplot::plot_grid(P1,P2,P3,ncol = 3,align = "vh")

#(3)beta----
D_before <- vegan::vegdist(map[,c(6:238)],method = 'bray')
pcoa5 <- cmdscale(D_before, k = 5, eig = TRUE)
eig=pcoa5$eig

set.seed(123)
ado = vegan::adonis(D_before ~ map$group,permutations = 999)
a = round(as.data.frame(ado$aov.tab[5])[1,1],3)
R2 = paste("adonis:R2 ",a, sep = "")
b = as.data.frame(ado$aov.tab[6])[1,1]
p_v = paste("p:",b, sep = "")
title1 = paste(R2,",",p_v, sep = "")
p0=ggplot(map, aes(x=V1, y=V2, color=group))  +#,shape=gender
  labs(x=paste0('PCoA'," 1 (",format(100*eig[1]/sum(eig),digits=4),"%)", sep=""),
       y=paste0('PCoA'," 2 (",format(100*eig[2]/sum(eig),digits=4),"%)", sep=""),
       title=title1) + geom_point(alpha=.9,size=2) + theme_classic() + 
  theme(text=element_text(face='bold', size=10),
        axis.text = element_text(face = 'bold',size = 10),
        plot.title=element_text(hjust=0.5),legend.position = 'none')+
  stat_ellipse(linetype=2,level=0.68,aes(group=group, colour=group))+
  scale_color_manual(values = group_col)+
  scale_y_continuous(limits = c(-0.5,0.5),breaks = seq(-0.5,0.5,0.25))+ 
  scale_x_continuous(limits = c(-0.5,0.5),breaks = seq(-0.5,0.5,0.25))+coord_fixed()

#(4)Enterotypes----
table(map$Enterotypes)
f=map[,c(6:238)]
colnames(f)=sapply(strsplit(colnames(f),'\\.'),function(x){x[6]})%>%gsub('g__','',.)

a1=colMeans(f[map$Enterotypes=='E_Ent',])[order(-colMeans(f[map$Enterotypes=='E_Ent',]))][1:8] 
#Escherichia
a2=colMeans(f[map$Enterotypes=='B_Ent',])[order(-colMeans(f[map$Enterotypes=='B_Ent',]))][1:8] 
#Bacteroides

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
temp = col_vector[5:27]
temp[c(1:7)]=col_vector[30:36]
col_vector=temp

mat=data.frame(genus=names(a1),relative=a1,Enterotype="E_Ent")
mat=rbind(mat,data.frame(genus=names(a2),relative=a2,Enterotype="B_Ent"))
mat$Enterotype=factor(mat$Enterotype,levels = c("E_Ent","B_Ent"))
mat$genus=factor(mat$genus,levels = unique(mat$genus))

P1=ggplot(mat,aes(x=Enterotype,y=relative,fill=genus))+
  geom_bar(stat="identity",position = "fill")+
  labs(x="",y="relative abundance (%)",title = 'Genus (top8)')+labs(fill="")+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5),text = element_text(size=10,face = 'bold'),
        axis.text = element_text(size=10,face = 'bold'))+#hjust = 0.5,vjust = 0.5
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values= col_vector[1:length(unique(mat$genus))])+
  coord_fixed(ratio = 3)

df1 = map %>% 
  count(group,Enterotypes) %>% 
  group_by(group) %>% 
  mutate(prop=prop.table(n))

P=ggplot(df1, aes(factor(group,levels =c('NC','SCS','SCD','MCI','AD')), prop*100, fill = Enterotypes))+
  geom_col() +
  geom_text(aes(label = scales::percent(prop)),
            position = position_stack(vjust = 0.5)) + 
  scale_fill_manual(values= col_vector[c(1,3)],limits=c('E_Ent','B_Ent'))+
  #scale_y_continuous(labels = scales::percent)+
  theme_classic()+
  theme(axis.text = element_text(face = 'bold',size = 10),legend.position = 'top',
        text= element_text(face = 'bold',size = 10))+labs(x='',y='Percentages (%)')

d=table(map$group,map$Enterotypes)
fisher.test(d)#Fisher’s exact test p-value = 0.1649

#3.Figure S2----
#(1)DA----
load(file = './01_taxa/all_marker.RData')
summary(results$g$log2FoldChange)
P1=ggplot(data=results$g,aes(y=gsub('g__','',taxa), x=group, fill=log2FoldChange)) + 
  geom_tile() + 
  scale_fill_gradientn(colours=c('#007A53','#009F4D',"#6CC24A", 'white',
                                 "#EFC06E", "#FFA300", '#BE5400'),
                       limits=c(-3.3,3.3))+theme_classic()+
  theme(panel.grid = element_blank(),text = element_text(face = 'bold',size = 9),
        axis.text = element_text(size=7,face='bold'),
        plot.title = element_text(hjust = 0.5,size = 9),
        axis.text.x=element_text(hjust = 1,vjust=0.5,angle=90)) +
  xlab('') + ylab('')+labs(title = 'LinDA: Genus-level')+coord_fixed()
#LinDA: P<0.05&abs(log2FoldChange)>1
