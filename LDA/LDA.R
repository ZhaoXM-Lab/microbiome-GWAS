# LDA.R 
rm(list=ls())

#1.data preparation----
map=read.csv('./01_taxa/feat_genus.csv',header=T)
otu=map[,c(6:238)]/100
rowSums(otu)
rownames(otu)=map$Row.names

library("OTUtable")
#abundance: The minimum threshold for percentage of reads attributed to a taxon in at least one sample
#persistence: The minimum threshold for the percentage of samples in which a taxon has been observed. 
t_filtered_abundance <- data.frame(filter_taxa(t(otu), abundance = 0.1, persistence = 10), check.names = F)
filtered_abundance <- data.frame(t(t_filtered_abundance))

filtered_otu <- filtered_abundance * 1000
filtered_otu <- sapply(round(filtered_otu), as.integer)
filtered_otu <- as.matrix(filtered_otu)

otu_names=data.frame(full_name=colnames(filtered_otu),genus=sapply(strsplit(colnames(filtered_otu),'\\.'),function(k){return(gsub('g__','',k[6]))}))
colnames(filtered_otu)=otu_names$genus
rownames(filtered_otu)=rownames(filtered_abundance)

#2.Hyperparameter selection + training----
library(slam)
library(topicmodels)
library(MetaTopics)

# set parameters 
dtm = as.simple_triplet_matrix(filtered_otu)
seed_num = 2022
fold_num = 5
kv_num = c(2, 5, 10, 15, 20, 25,30,35,40) #topic numbers to compare 
sp = smp(cross = fold_num, n = nrow(dtm), seed = seed_num) 
control = list(seed = seed_num, burnin = 100, thin = 10, iter = 100) 

ctmK = selectK(dtm = dtm, kv = kv_num, SEED = seed_num,
               cross = fold_num, sp = sp, method = "Gibbs",
               control = control)

# compare results
plot_perplexity(ctmK, kv_num)
dev.off()
# log―likihood+perplexity
# save(ctmK,file='./LDA_model_K.RData')

load('./LDA_model_K.RData')
n_topic=15
seed_num = 666
Gibbs_model = LDA(dtm, k = n_topic, method = "Gibbs", control = list(seed = seed_num, burnin = 1000, thin = 100, iter = 1000))
# save(Gibbs_model,file='./Gibbs_model.RData')

