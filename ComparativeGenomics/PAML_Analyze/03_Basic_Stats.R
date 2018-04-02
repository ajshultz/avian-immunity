setwd("~/Dropbox/BirdImmuneGeneEvolution")
library(tidyverse)

######################################################################
##Basic dataset characteristics#######################################
######################################################################
#How many genes in gene tree and species tree datasets?
all_res_gene_ncbi %>%
  summarize(n())
all_res_sp_ncbi %>%
  summarize(n())


#How many hogs could be assigned a chicken gene ID?
all_res_gene_ncbi %>%
  filter(!is.na(entrezgene)) %>%
  summarize(n())

#How many hogs could be assigned a zebra finch gene ID?
all_res_gene_ncbi %>%
  filter(!is.na(entrezgene_zf)) %>%
  summarize(n())

#How many hogs could be assigned to both a chicken or zebra finch gene ID?
all_res_gene_ncbi %>%
  filter(!is.na(entrezgene), !is.na(entrezgene_zf)) %>%
  summarize(n())

#How many hogs coudl not be assigned to a gene ID?
all_res_gene_ncbi %>%
  filter(is.na(entrezgene), is.na(entrezgene_zf)) %>%
  summarize(n())



all_res_gene_ncbi <- all_res_gene_ncbi %>%
  filter(!is.na(pval_busted) & !is.na(PVal_m1m2) & !is.na(PVal_m2m2a) & !is.na(PVal_m7m8) & !is.na(PVal_m8m8a) & !is.na(total_sel.n))

all_res_sp_ncbi <- all_res_sp_ncbi %>%
  filter(!is.na(pval_busted) & !is.na(PVal_m1m2) & !is.na(PVal_m2m2a) & !is.na(PVal_m7m8) & !is.na(PVal_m8m8a) & !is.na(total_sel.n))

save(all_res_gene_ncbi,all_res_sp_ncbi,file="02_output_annotated_data/all_res_ncbi.Rdat")




#How many genes in gene tree and species tree datasets?
all_res_gene_ncbi %>%
  summarize(n())
all_res_sp_ncbi %>%
  summarize(n())

#either chicken or zebra finch?
all_res_gene_ncbi %>%
  filter(!is.na(entrezgene_zf) | !is.na(entrezgene)) %>%
  summarize(n())

all_res_gene_zf_hs %>%
  filter(!is.na(entrezgene)) %>%
  filter(!is.na(ensembl_gene_id))




#############################################################################
#############Calculate some general stats on the dataset#####################
#############################################################################

stats_table <- matrix(nrow=2, ncol=9)

stats_table[,1] <- c("gene","species")

#How many genes with all tests?
stats_table[1,2] <- all_res_gene_ncbi %>%
  filter(!is.na(pval_busted) & !is.na(PVal_m1m2) & !is.na(PVal_m2m2a) & !is.na(PVal_m7m8) & !is.na(PVal_m8m8a) & !is.na(total_sel.n)) %>%
  summarise(n()) %>% pull
stats_table[2,2] <- all_res_sp_ncbi %>%
  filter(!is.na(pval_busted) & !is.na(PVal_m1m2) & !is.na(PVal_m2m2a) & !is.na(PVal_m7m8) & !is.na(PVal_m8m8a) & !is.na(total_sel.n)) %>%
  summarise(n()) %>% pull

#How many hogs are selected according to each test?


#How many hogs are selected in all PAML tests
all_res_gene_ncbi %>%
  filter(FDRPval_m2m2a < 0.05, FDRPval_m8m8a < 0.05, FDRPval_m1m2 < 0.05, FDRPval_m7m8 < 0.05)
all_res_sp_ncbi %>%
  filter(FDRPval_m2m2a < 0.05, FDRPval_m8m8a < 0.05, FDRPval_m1m2 < 0.05, FDRPval_m7m8 < 0.05)

all_selected_hogs <- hyphy_paml_res[hyphy_paml_res[,"FDRPval_m2m2a"]<0.05 & hyphy_paml_res[,"FDRPval_m8m8a"]<0.05 & hyphy_paml_res[,"FDRPval_m1m2"]<0.05 & hyphy_paml_res[,"FDRPval_m7m8"]<0.05 & hyphy_paml_res[,"FDRPval_busted"]<0.05,"hog"]
all_selected_hogs <- as.character(all_selected_hogs[!is.na(all_selected_hogs)])
length(all_selected_hogs)

#How many branches are selected by BS-Rel for HOGs selected by all tests
hist(hyphy_paml_res[all_selected_hogs,"total_sel.s"],breaks=20)
hist(hyphy_paml_res[all_selected_hogs,"total_sel.n"],breaks=20)
summary(hyphy_paml_res[all_selected_hogs,"total_sel.s"])
summary(hyphy_paml_res[all_selected_hogs,"total_sel.n"])