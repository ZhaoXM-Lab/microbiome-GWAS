#02_ES-Ana_associated_with_AD.R 

rm(list=ls())

library(topicmodels)
library(MetaTopics)
library(dplyr)

load('./Gibbs_model.RData')

#1.ES-Ana probability: logistic----
map=read.csv('./03_table/00_meta_all.csv',header=T)
map$age=map$Age

load('./Gibbs_model.RData')
n_topic=15
LDA_20_1_df <- as.data.frame(Gibbs_model@gamma)
names(LDA_20_1_df) <- paste0('Topic_',1:n_topic)
LDA_20_1_df$sample <- Gibbs_model@documents
Nutribiome_Topics <- merge(map[,c('FID','group','age','gender','MMSE','MoCA_B')], LDA_20_1_df, 
                           by.x = "FID",by.y = 'sample', sort = F)

lg.models=list()
results=list()
data.sig.confusion=list()

for(i in c('AD')){
  train.data=subset(Nutribiome_Topics[,c('group',paste0('Topic_',1:15))],group%in%c(i,'NC'))%>%
    select(group,Topic_8)
  train.data$group=ifelse(train.data$group=='NC',0,1)%>%factor(levels = c(0,1))
  lg.models[[i]] <- glm(group~., data=train.data, family='binomial')
  results[[i]]=pROC::roc(response=train.data$group,predictor=lg.models[[i]]$fitted.values,
                         levels = c(0,1),
                         direction='<')
  predictions.class <- ifelse(lg.models[[i]]$fitted.values>=0.5,1,0)%>%factor(levels = c(0,1))
  data.sig.confusion[[i]] <- caret::confusionMatrix(predictions.class, train.data$group)
  train.data=NULL
}
results[['AD']]$auc
pROC::ci(results[['AD']])
data.sig.confusion[['AD']]

#ROC curve----
#Discovery cohort: eval1
#Japan cohort: eval2
#Replication cohort: eval3

plot(
  NULL,
  xlim = c(0, 1),
  ylim = c(0, 1),
  xlab = "False positive rate",
  ylab = "True positive rate",
  type = "n"
)
title(paste("ROC curve", sep = " "))
abline(a = 0, b = 1, lty = 3)

lines(1 - eval1$specificities, eval1$sensitivities,
      col = 'red', lwd = 2)
lines(1 - eval2$specificities, eval2$sensitivities,
      col = 'orange', lwd = 2)
lines(1 - eval3$specificities, eval3$sensitivities,
      col = 'skyblue', lwd = 2)

# plot CI
x = as.numeric(rownames(eval1$ci))
yl = eval1$ci[,1]
yu = eval1$ci[,3]
polygon(1 - c(x, rev(x)), c(yl, rev(yu)),
        col = alpha('red', alpha=0.1),
        border = NA)

x = as.numeric(rownames(eval2$ci))
yl = eval2$ci[, 1]
yu = eval2$ci[, 3]
polygon(1 - c(x, rev(x)), c(yl, rev(yu)),
        col = alpha('orange', alpha=0.2),
        border = NA)

x = as.numeric(rownames(eval3$ci))
yl = eval3$ci[, 1]
yu = eval3$ci[, 3]
polygon(1 - c(x, rev(x)), c(yl, rev(yu)),
        col = alpha('skyblue', alpha=0.3),
        border = NA)

legend("bottomright",legend =c(paste0('AUC (train on Discovery): ',round(eval1$auc,digits = 2)),
                               paste0('AUC (test on Japan): ',round(eval2$auc,digits = 2)),
                               paste0('AUC (test on Replication): ',round(eval3$auc,digits = 2))),
       col=c('red','orange','skyblue'),lty = c(1,1,1),cex=0.8)
dev.off()

#2.key genera of ES-Ana----
#2.1 fold change----
my_beta <- Gibbs_model@beta
colnames(my_beta) <- Gibbs_model@terms
rownames(my_beta) <- paste("Topic", 1:nrow(my_beta),sep = "")
my_beta <- exp(my_beta)

xx=sort(my_beta[8,],decreasing = T)
xx=names(xx[xx>=0.01])#8 genera

meta=read.csv(file = './03_table/00_meta_all.csv',header=T)
rownames(meta)=meta$FID
map=read.csv('./01_taxa/feat_genus.csv',header=T)
stopifnot(meta$FID==map$Row.names)

otu=map[,grepl('k__',colnames(map))]/100#233
rownames(otu)=map$Row.names
otu_names=data.frame(full_name=colnames(otu),
                     genus=sapply(strsplit(colnames(otu),'\\.'),function(k){return(gsub('g__','',k[6]))}))
colnames(otu)=otu_names$genus
otu=otu[,xx]#8

my_FC=function(f,m,case,control,group){
  stopifnot(rownames(f)==rownames(m))
  mean_1=colMeans(f[m[,group]==case,])#group1:case
  mean_2=colMeans(f[m[,group]==control,])#group2:control
  log2fc=log2(mean_1/mean_2)
  
  return(log2fc)
}

results=data.frame(taxa=colnames(otu),
                   FC.AD=my_FC(otu,meta,'AD','NC','group'),
                   FC.MCI=my_FC(otu,meta,'MCI','NC','group'),
                   FC.SCD=my_FC(otu,meta,'SCD','NC','group'),
                   FC.SCS=my_FC(otu,meta,'SCS','NC','group'))
otu_names=subset(otu_names,genus%in%colnames(otu))
rownames(results)=results$taxa
results=results[,c('FC.SCS','FC.SCD','FC.MCI','FC.AD')]

library(ComplexHeatmap)
col_fun = circlize::colorRamp2(c(-2,-1,-0.5,0,0.5,1,2),
                               c('#007A53','#009F4D',"#6CC24A", 'white',
                                 "#EFC06E", "#FFA300", '#BE5400'))


otu_names$phylum=sapply(strsplit(otu_names$full_name,'\\.'),function(k){return(gsub('p__','',k[2]))})
otu_names$family=sapply(strsplit(otu_names$full_name,'\\.'),function(k){return(gsub('f__','',k[5]))})
rownames(otu_names)=otu_names$genus
otu_names=otu_names[colnames(otu),]

left_anno=otu_names[,c('phylum','family')]
left_annos = HeatmapAnnotation(df = left_anno,which = "row",
                               col=list(
                                 phylum=c(Firmicutes="#E6DCEC",Actinobacteria="#CEE4EF"),
                                 family=c(Lachnospiraceae='#CAB2D6',Ruminococcaceae='#F9E9A4',Eggerthellaceae='#96C1D7')
                               ))

Heatmap(results,col=col_fun, width = unit(4, "cm"), height = unit(8, "cm"),
        column_title = 'Topic 8',
        name = "log2(FC) (VS NC)", 
        row_split = 3,cluster_rows = T,cluster_columns=F,
        left_annotation=left_annos,
)

#2.2 importance rank in the rf model----
#define functions:----
#(1)my_rf_test
#(2)my_im_plot

my_rf_test=function(f,m,group,group_level,
                    p0=NULL,seed=123,do.rfe=T,marker=NULL,subsets=c(10,15,20,30,40,50),
                    times=NULL,k=NULL,level='s',str_t='\\|',cut_name=F,
                    norm.param=NULL,type='taxa',covar=T,raw=T,do_smote=F){
  require(caret)
  library(randomForest)
  
  library(magrittr)
  library(dplyr)
  library(purrr)
  
  colnames(m)[colnames(m)==group]='Group'
  group='Group'

  meta=merge(f,m,by='row.names',sort = F)
  rownames(meta)=meta$Row.names
  meta=meta[,c("gender","age","Group")]
  
  meta=subset(meta,Group%in%group_level)
  f=f[rownames(meta),]
  f=f[,colSums(f)>0]
  
  #1.data preparations
  if(raw==T){
    f.filtered=f[,grepl('p__|c__|o__|f__|g__|s__',colnames(f))]/100
  }else{
    #TSS: genus
    f.filtered=f
  }
  
  #2.normalization: log.std by default
  my_normal=function(f,method,log.n0 = 1e-06,
                     sd.min.q = 0.1,norm.param=NULL){
    
    par=NULL
    if(method=='TSS'){
      f=sweep(f,1,rowSums(f),`/`)%>%as.data.frame()#row:sample
    }
    
    rowSds=function(x){
      apply(x, 1, sd)
    }
    if (method == "log.std") {
      if(is.null(norm.param)){
        f=t(f)#row:genus
        feat.log <- log10(f + log.n0)
        m <- rowMeans(feat.log)
        s <- rowSds(feat.log)
        q <- quantile(s, sd.min.q, names = FALSE)
        stopifnot(q > 0)
        f <- (feat.log - m) / (s + q)
        f=t(f)
        par=list()
        par$feat.mean <- m
        par$feat.adj.sd <- s + q
      }else{
        feat.log <- log10(f + log.n0)
        feat.norm <- (feat.log - norm.param$feat.mean) /
          norm.param$feat.adj.sd
        f=t(feat.norm)
      }
    }
    
    if(method == "log"){
      f=log10(f + 1)
    }
    
    return(list(f=f,norm.param=par))
  }
  
  norm.train=my_normal(f.filtered,method = 'log.std',norm.param=norm.param)
  f.norm=norm.train[['f']]
  
  #3.sample spliting
  if (cut_name==T) {
    if(level=='s'){
      colnames(f.norm)=sapply(strsplit(colnames(f.norm),str_t),function(x){x[[7]]})%>%gsub('s__','',.)
    }else{
      colnames(f.norm)=sapply(strsplit(colnames(f.norm),str_t),function(x){x[[6]]})%>%gsub('g__','',.)
    }
  }
  stopifnot(rownames(f.norm)==rownames(f))
  
  if(type=='taxa'){
    f.norm=f.norm
  }else if(type=='SNP'){
    f.norm=f[,!grepl('p__|c__|o__|f__|g__|s__',colnames(f))]
  }else{
    #both
    f[,grepl('p__|c__|o__|f__|g__|s__',colnames(f))]=f.norm
    f.norm=f
  }
  
  f.norm=as.data.frame(f.norm)
  
  #4.confounder: age+sex
  stopifnot(rownames(f)==rownames(meta))
  if(covar==T){
    f.norm$age=meta$age
    f.norm$gender=ifelse(meta$gender=='male',1,0)
  }
  
  feature_name=data.frame(model=paste0("f",1:ncol(f.norm)),full=colnames(f.norm))
  colnames(f.norm)=feature_name$model
  
  m=meta
  stopifnot(rownames(f.norm)==rownames(m))
  f.norm=as.data.frame(f.norm)
  
  dat=cbind(f.norm,Group=m[,group])%>%as.data.frame()
  dat[,group]=factor(dat[,group],levels = group_level)
  
  #5.cross validation
  set.seed(seed)
  sample.split <- createMultiFolds(dat[,group], times = times, k=k)
  rf.models <- list()
  pred.matrix <- list()
  confusion.result <- list()
  
  #5.1
  for (iter in names(sample.split)) {
    train.data  <- dat[sample.split[[iter]], ]
    test.data <- dat[-sample.split[[iter]], ]
    set.seed(12345)
    rf.models[[iter]] <- randomForest(Group ~ .,train.data, proximity = T, importance = T, ntree = 500)

    predictions <- predict(rf.models[[iter]], test.data, type = "prob")
    pred.matrix[[iter]] <- predictions
    
    predictions.class <- predict(rf.models[[iter]], test.data, type= "response")
    confusion.result[[iter]] <- predictions.class
  }
  
  #5.2
  pred.matrix.list <- list()
  pred.df <- data.frame()
  pred.df.2 <- data.frame()
  
  fold.str <- paste("Fold", gsub(" ", "0", format(1:k)), sep = "")
  rep.str <- paste("Rep", gsub(" ", "0", format(1:times)), sep = "")
  for (rep in rep.str) {
    for (fold in fold.str) {
      iter <- paste(fold, rep, sep = ".")
      pred.df.2 <- data.frame(pred.matrix[[iter]], stringsAsFactors = F)
      pred.df <- rbind(pred.df, pred.df.2)
    }
    for (grp in names(pred.df)) {
      
      y <- data.frame(pred.df[,grp], stringsAsFactors = F);
      names(y) <- grp;
      y$Sample <- rownames(pred.df)
      pred.matrix.list[[grp]][[rep]] <- y
    }
    pred.df <- data.frame();
  }
  
  ## --- importance scores --- ##
  importance <- lapply(rf.models, function(data){ imp <- importance(data, type=1)})
  importance.df <- purrr::reduce(importance, cbind)
  importance.df <- data.frame(importance.df, stringsAsFactors = F)
  names(importance.df) <- names(rf.models)
  data=list(sample.split = sample.split,
            rf.models = rf.models,
            pred.matrix = pred.matrix,
            confusion.result = confusion.result,
            importance.df = importance.df,
            pred.matrix.list = pred.matrix.list
  )
  
  return(list(data=data,feature_name=feature_name))
}

my_im_plot=function(importance.df,top=20,title,
                    f=NULL,feature_name=NULL){
  library(dplyr)
  library(ggplot2)
  if(!is.null(feature_name)){
    stopifnot(rownames(importance.df)==feature_name$model)
    rownames(importance.df)=feature_name$full
  }
  
  feat.relaimpo <- importance.df %>% 
    rowMeans(.) %>% 
    tibble( relaimpo = ., species = names(.) ) %>% 
    arrange( desc (relaimpo ) );
  feat.select <- feat.relaimpo %>% dplyr::slice( 1:top )
  
  
  ## -- 2. relative importance plot --
  #k*n*top
  feat.importance.select <- as.data.frame( t( importance.df[ feat.select$species, ] ) ) %>% tidyr::gather( feature, relaimpo );
  
  feat.importance.select.plot <- 
    feat.importance.select %>% ggplot( aes( reorder( feature, relaimpo, median) ,  relaimpo) ) + 
    geom_boxplot( aes( fill = feature ) ) + coord_flip() + theme_classic()+
    theme( legend.position = "none" ,
           axis.text = element_text(face = 'bold',size = 10),
           text= element_text(face = 'bold',size = 10),
           plot.title = element_text(hjust=0.5))+
    labs(y='importance',x='feature',title=title)
  
  return(list(f.plot=feat.importance.select.plot,f.sel=feat.select$species))
}

#define functions: over----

load('./03_table/all_WGS_16S_data_B_F.RData')
model0=my_rf_test(f=datas$Bacteria$our1$feat,m=datas$Bacteria$our1$meta,
                  group='group',group_level=c('NC','AD'),
                  norm.param=NULL,seed=123,times=3,k=5,type='taxa',covar=F,raw=F)

ca_p=my_im_plot(importance.df=model0$data$importance.df,top=20,
                title='AD-NC',feature_name = model0$feature_name)
ca_p$f.plot
