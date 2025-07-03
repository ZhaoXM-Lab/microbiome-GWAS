# Topic_distribution.R 

rm(list=ls())
library(topicmodels)
library(MetaTopics)
library(ggplot2)

load('./Gibbs_model.RData')
map=read.csv('./01_taxa/feat_genus.csv',header=T)

#1.topic-genus probability----
P1=plot_beta(Gibbs_model,prob=0.01)

#pie chart: Topic_8----
my_beta <- Gibbs_model@beta
colnames(my_beta) <- Gibbs_model@terms
rownames(my_beta) <- paste("Topic", 1:nrow(my_beta),sep = "")
my_beta <- exp(my_beta)

xx=sort(my_beta[8,],decreasing = T)
sum(xx[xx>=0.01])#0.9602766
xx[xx>=0.01]
sum(xx)
xx=xx[xx>=0.01]
xx=c(xx,1-sum(xx))
names(xx)[9]='Others'
xx=data.frame(genus=names(xx),percent=xx)

xx$percent=round(xx$percent*100)
cbPalette <- c("#6A3D9A","#936FB3","#8891DB","#CC88B0","#D2ACB6","#1F78B4","#7B97B1","#DAC256",'#B5B5B6')
names(cbPalette)=xx$genus

xx$labs <- paste(xx$genus, ' (', xx$percent, '%)', sep = '')
xx$genus=factor(xx$genus,levels=xx$genus)

library(ggpubr)
p3 <- ggpie(xx, 'percent',  #绘图，只用写频数就行，切记不用再写分组
            fill = 'genus', label = 'labs',
            palette = cbPalette,
            lab.font = c(4),lab.pos = 'out')+labs(title = 'Genera probability distribution in Topic_8')+
  theme(legend.position = "none") #设置标签，标签的大小为4， 
p3

#2.sample-topic probability-DM----
#define function: ----
#(1)my_DM()
#(2)matrix_barplot

my_DM=function(Nutribiome_Topics=NULL,n_topic,
               sd.var=NULL,Group=NULL,Group.level=NULL,sep='-'){
  library(DirichletReg)
  library(dplyr)
  
  if(is.null(sep)){
    Topics_matrix <- DR_data(Nutribiome_Topics[, paste0('Topic',1:n_topic)], trafo = F)
  }else{
    Topics_matrix <- DR_data(Nutribiome_Topics[, paste0('Topic_',1:n_topic)], trafo = F)
  }
  
  Nutribiome_Topics$gender=relevel(factor(Nutribiome_Topics$gender),ref='female')

  (sd.age <- sd(Nutribiome_Topics$age, na.rm = T)) 
  Nutribiome_Topics$age.sd <- Nutribiome_Topics$age/sd.age
  
  if(!is.null(sd.var)){
    Nutribiome_Topics$var=Nutribiome_Topics[,sd.var]/sd(Nutribiome_Topics[,sd.var], na.rm = T)
    model <- DirichReg(Topics_matrix ~ var + age.sd + gender,
                       data = Nutribiome_Topics)
    
    model.sum <- (summary(model)) 
    coef2 <- model.sum$coef.mat 
    selected <- grep("var", dimnames(coef2)[[1]]) 
    coefmatrix <- coef2[selected, ]

    estvector <- coefmatrix[, "Estimate"]
    estmatrix <- matrix(estvector, ncol = n_topic,byrow = F)
    if(is.null(sep)){
      colnames(estmatrix) <-paste0('Topic',1:n_topic) 
    }else{
      colnames(estmatrix) <-paste0('Topic_',1:n_topic) 
    }
    rownames(estmatrix) = sd.var
    
    pvector <- coefmatrix[, "Pr(>|z|)"]
    pmatrix <- matrix(pvector, ncol = n_topic,byrow = F)
    
    if(is.null(sep)){
      colnames(pmatrix) <-paste0('Topic',1:n_topic) 
    }else{
      colnames(pmatrix) <-paste0('Topic_',1:n_topic) 
    }
    rownames(pmatrix) <- sd.var
  }else{

    stopifnot(!is.null(Group))

    Nutribiome_Topics$Group=factor(Nutribiome_Topics[,Group],levels = Group.level)
    model <- DirichReg(Topics_matrix ~ Group + age.sd + gender,
                       data = Nutribiome_Topics)
    
    model.sum <- (summary(model)) 
    coef2 <- model.sum$coef.mat 
    selected <- grep("Group", dimnames(coef2)[[1]]) 
    coefmatrix <- coef2[selected, ]

    estvector <- coefmatrix[, "Estimate"]
    estmatrix <- matrix(estvector, ncol = n_topic,byrow = F)
    if(is.null(sep)){
      colnames(estmatrix) <-paste0('Topic',1:n_topic) 
    }else{
      colnames(estmatrix) <-paste0('Topic_',1:n_topic) 
    }

    rownames(estmatrix) <- Group.level[-1]
    
    pvector <- coefmatrix[, "Pr(>|z|)"]
    pmatrix <- matrix(pvector, ncol = n_topic,byrow = F)
    if(is.null(sep)){
      colnames(pmatrix) <-paste0('Topic',1:n_topic) 
    }else{
      colnames(pmatrix) <-paste0('Topic_',1:n_topic) 
    }
    rownames(pmatrix) <- Group.level[-1]
  }
  return(list(model=model,p=pmatrix,beta=estmatrix))
}

matrix_barplot = function(data, group_by=NULL, pvals=NULL, xlab='', ylab='Frequency', 
                          value='mean', error='se', legend.title='Groups', 
                          colors='Paired', pos='dodge', border=NA,
                          coord_flip=FALSE, sig_only=F, do.facet=F,
                          Group.level =c('NC','SCS','SCD','MCI','AD')){
  library(data.table)
  library(tidyr)
  library(ggplot2)


  if(is.null(group_by)){group_by = rownames(data)}
  if(nlevels(group_by) == 0){group_by = factor(group_by,levels = Group.level)}
  
  # Select significant comparisons
  if(sig_only == TRUE){
    j = apply(pvals, 2, min) <= .05
    if(sum(j) == 0){return(NULL)}
    data = data[,j,drop=F]
    pvals = pvals[,j,drop=F]
  }
  
  # Construct input data
  names = colnames(data)#topic
  data = data.frame(group=group_by, data)
  group_levels = levels(group_by)
  
  data = tidyr::gather(data,key = 'x',value = 'y',-c(group))#252*20
  
  # Value function
  if(value == 'mean'){vf = mean} else if(value == 'median'){vf = median} else {stop()}
  
  # Error function
  se = function(x, na.rm=T){sd(x, na.rm=na.rm)/sqrt(length(x))}
  if(error == 'sd'){ef = sd} else if(error == 'se'){ef = se} else {ef = function(x, ...){0}}
  
  # Estimate error bars
  data=data%>%group_by(group,x)%>%mutate(u=vf(y, na.rm=T),s=ef(y, na.rm=T))%>%ungroup()
  data=as.data.table(data)
  
  # Add p-values 
  if(!is.null(pvals)){
    pvals = as.data.frame(pvals) %>% tibble::rownames_to_column('group') %>% gather(x, pval, -group)
    pvals = as.data.table(pvals)
    
    setkeyv(pvals,c('group','x'))
    setkeyv(data,c('group','x'))
    
    data = merge(data, pvals, by=c('group','x'),all =T)
    data$lab1=dplyr::case_when(data$pval <=0.001 ~ '***',
                               data$pval <=0.01 ~ '**',
                               data$pval <=0.05 ~ '*',
                               data$pval <=0.1 ~ '.',
                               data$pval >0.1 ~ '')
    
  }
  
  if(coord_flip == TRUE){names = rev(names); group_levels=rev(group_levels)}
  data$x = factor(data$x, levels=names)    
  data$group = factor(data$group, levels=group_levels)
  
  # Get colors
  if(length(colors) == 1){colors = set.colors[1:length(group_levels)]}
  
  # Plot data
  if(pos == 'stack'){
    p = ggplot(data) + geom_bar(aes(x=x, y=u, fill=group), colour=border, size=.25, stat='identity')
    if(error %in% c('sd', 'se')){p = p + geom_errorbar(aes(x=x, ymin=u-s, ymax=u+s, fill=group), stat='identity', width=.25)}
  } else {
    pos = position_dodge(.9)
    p = ggplot(data) + geom_bar(aes(x=x, y=u, fill=group), colour=border, size=.25, stat='identity', position=pos)
    if(error %in% c('sd', 'se')){
      p = p + geom_errorbar(aes(x=x, ymin=u-s, ymax=u+s,fill=group), stat='identity', position=pos, width=.25)
    }
  }
  
  p = p + 
    scale_fill_manual(values=colors, name=legend.title) + xlab(xlab) + ylab(ylab) +
    scale_color_manual('', values=c('#000000', '#999999', '#cccccc'), guide='none')
  
  # Facet wrap
  if(do.facet == TRUE){
    p = p + facet_grid(group ~ ., scales='free')
  }
  
  dy = max(data$u + data$s, na.rm=T)*.01
  if(coord_flip == FALSE){
    p = p + theme_classic()+theme(axis.text.x = element_text(angle = -45, hjust = 0),text = element_text(face = 'bold'))
    if(!is.null(pvals)){
      p = p + geom_text(aes(x=x, y=u+s+dy, label=lab1, group=group), hjust='center', vjust=0, size=5, angle=0, position=pos)
    }
  } else {
    p = p + coord_flip()
    if(!is.null(pvals)){
      p = p + geom_text(aes(x=x, y=u+s+dy, label=lab1, group=group), hjust='center', vjust=1, size=5, angle=90, position=pos)
    }
  }
  return(p)
}


#define function: over----

n_topic=15
library(DirichletReg)
LDA_20_1_df <- as.data.frame(Gibbs_model@gamma)
names(LDA_20_1_df) <- paste0('Topic_',1:n_topic)
LDA_20_1_df$sample <- Gibbs_model@documents
Nutribiome_Topics <- merge(map[,c(1:5)], LDA_20_1_df, by.x = "Row.names",by.y = 'sample', sort = F)
ep.pct = 100*Nutribiome_Topics[, paste0('Topic_',1:n_topic)]/rowSums(Nutribiome_Topics[, paste0('Topic_',1:n_topic)])

D1=my_DM(Nutribiome_Topics,15,Group = c('group'),Group.level =c('NC','SCS','SCD','MCI','AD'))
Nutribiome_Topics$group=factor(Nutribiome_Topics$group,levels = c('NC','SCS','SCD','MCI','AD'))

group_col=c("#982b2b","#db6968","#EDB3B3","#D0E5D0","#459943")
names(group_col)=c('AD','MCI','SCD','SCS','NC')

#2.1 AD vs. NC----
D1$p#0.04192971
ids=Nutribiome_Topics$group%in%c('AD','NC')
pmatrix=D1$p['AD',]%>%t()%>%as.matrix()
rownames(pmatrix)='AD'
p1 = matrix_barplot(data=ep.pct[ids,], group_by=Nutribiome_Topics$group[ids], pvals=pmatrix, colors=group_col)

#2.2 all diagnosis groups----
rowSums(Nutribiome_Topics[, paste0('Topic_',1:n_topic)])
summary(Nutribiome_Topics[, paste0('Topic_',1:n_topic)])

Nutribiome_Topics[, paste0('Topic_',1:n_topic)]=Nutribiome_Topics[, paste0('Topic_',1:n_topic)]*100
data = reshape2::melt(Nutribiome_Topics[, c('group',paste0('Topic_',1:n_topic))],id=c('group'),variable.name='topic',value.name = 'frequency')
map=data%>%group_by(group,topic)%>%summarise(median=median(frequency))%>%ungroup()%>%as.data.frame()
map$lower_bound=NA
map$upper_bound=NA

# bootstrapping
library(boot)
median_fun <- function(data, indices) {
  median(data[indices])
}

for( i in 1:nrow(map)){
  set.seed(123)
  bootstrap_result <- boot(data=data$frequency[data$group==map$group[i]&data$topic==map$topic[i]], statistic=median_fun, R=1000)
  xx=boot.ci(bootstrap_result, type="perc", conf=0.95)# 95% CI
  map[i,c('lower_bound','upper_bound')]=xx$percent[1,c(4,5)]
}

map$group=factor(map$group,levels = c('NC','SCS','SCD','MCI','AD'))
P1=ggplot(map, aes(y = median, x = group, fill = group)) + 
  facet_grid(~topic, scales="free", space="free")+
  geom_errorbar(aes(x = group, ymin=lower_bound, ymax=upper_bound, color=group), stat='identity', width=.25) + 
  geom_point(shape = 21, size = 3, stroke = 0) + 
  stat_summary(fun = "median", geom = "line", aes(group = topic))+ 
  theme_classic() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+
  labs(x = "", y = "Frequency (%)") + 
  scale_fill_manual(values=group_col)+
  scale_color_manual(values=group_col)

#2.3 Dirichlet regression----
samples=read.csv('./01_taxa/11_Topic_sample_v2.csv',header = T)
Topics <- merge(samples[,c('FID','MMSE','MoCA_B')], Nutribiome_Topics, by.x = "FID",by.y = 'Row.names', sort = F)
D4=my_DM(Topics,15,sd.var='MMSE')
D5=my_DM(Topics,15,sd.var='MoCA_B')
D4$p#0.04696557
D4$beta#0.1133527
D5$p#0.04267945
D5$beta#0.1174418

beta=rbind(D1$beta,D4$beta,D5$beta)
max(abs(beta))
p=rbind(D1$p,D4$p,D5$p)

library(corrplot)
col <- colorRampPalette(c("royalblue4", "royalblue2", "royalblue1", "White", "orangered1", "orangered2", "orangered4"), interpolate = "linear")

max(beta)
min(beta)

corrplot(beta, is.corr = FALSE, method = c("circle"), addgrid.col = NA, 
         col = col(200), tl.cex = 0.6, tl.col = "black",cl.pos = "r",
         cl.length = 9,col.lim = c(-0.4,0.4),
         p.mat=p, sig.level = 0.05, insig = "label_sig")

