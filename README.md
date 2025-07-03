# Introduction

The code is used for the work of “Multi-omics analysis reveals associations between gut microbial subgroup alterations and host genetic aberrations in Alzheimer’s disease ”. This code includes two parts, the first part is the LDA topic model, and the second part is about statistical analyses.



# LDA

* LDA.R, the main process for identifying microbial subgroups by LDA topic model 
* Topic_distribution.R, the topic-sample probability and the genus-topic probability matrix of the fitted LDA model; Cognitive score-topic and disease-topic associations



# Analysis

* 01_statistic_of_gut_microbiome.R, visualizations of gut metagenomic characteristics including alpha-/beta-diversity, enterotype distribution, and differential abundance genera.
* 02_ES-Ana_associated_with_AD.R, the importance of ES-Ana probability and 8 key genera in distinguishing AD with NC
* 03_GWAS_for_ES-Ana.R, GWAS results for ES-Ana probability in the Discovery cohort and meta-analysis 
* 04_PRS_association.R, associations between PRS of AD and ES-Ana probability or cognitive scores
* 05_mediation.R, mediation analysis inferring the putative SNP-microbe-AD routes
* 06_enrichment_analysis.R, pheWAS results and enriched biological functions and traits



# Contact information

If you have any questions about this code, please contact the corresponding author (xmzhao@fudan.edu.cn) or the first author (21110850029@m.fudan.edu.cn).

