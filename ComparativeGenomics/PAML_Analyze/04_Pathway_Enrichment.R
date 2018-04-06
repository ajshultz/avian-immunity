library(tidyverse)
library(DOSE)
library(clusterProfiler)
library(pathview)

setwd("~/Dropbox/BirdImmuneGeneEvolution/")

#######################################################################################################################
#Get significant and tested gene lists for chicken
#######################################################################################################################

#load processed, NCBI annotated data
load("02_output_annotated_data/all_res_ncbi.Rdat")

#Get list of entrezgenes significant in all tests, and lists of genes present in all tests
all_tested_gene <- all_res_gene_ncbi %>%
  filter(!is.na(FDRPval_m1m2) & !is.na(FDRPval_m2m2a) & !is.na(FDRPval_m7m8) & !is.na(FDRPval_m8m8a) & !is.na(FDRPval_busted)) %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)

all_sig_gene <- all_res_gene_ncbi %>%
  filter(!is.na(FDRPval_m1m2) & !is.na(FDRPval_m2m2a) & !is.na(FDRPval_m7m8) & !is.na(FDRPval_m8m8a) & !is.na(FDRPval_busted)) %>%
  filter(!is.na(entrezgene)) %>%
  filter(FDRPval_m1m2 < 0.05 & FDRPval_m2m2a < 0.05 & FDRPval_m7m8 < 0.05 & FDRPval_m8m8a < 0.05 & FDRPval_busted < 0.05) %>%
  pull(entrezgene)

all_tested_sp <- all_res_sp_ncbi %>%
  filter(!is.na(FDRPval_m1m2) & !is.na(FDRPval_m2m2a) & !is.na(FDRPval_m7m8) & !is.na(FDRPval_m8m8a) & !is.na(FDRPval_busted)) %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)

all_sig_sp <- all_res_sp_ncbi %>%
  filter(!is.na(FDRPval_m1m2) & !is.na(FDRPval_m2m2a) & !is.na(FDRPval_m7m8) & !is.na(FDRPval_m8m8a) & !is.na(FDRPval_busted)) %>%
  filter(!is.na(entrezgene)) %>%
  filter(FDRPval_m1m2 < 0.05 & FDRPval_m2m2a < 0.05 & FDRPval_m7m8 < 0.05 & FDRPval_m8m8a < 0.05 & FDRPval_busted < 0.05) %>%
  pull(entrezgene)

#How many genes are different between the two significant sets?
length(setdiff(all_sig_sp,all_sig_gene))

#######################################################################################################################
#Pathway Enrichment
#######################################################################################################################

#Chicken pathway enrichemnt tests:
#Gene trees, p<0.05, q<0.05
all_genes_k <- enrichKEGG(all_sig_gene,organism="gga",pvalueCutoff=0.05,pAdjustMethod="BH",qvalueCutoff=0.05,universe=all_tested_gene,keyType="ncbi-geneid")

write_csv(summary(all_genes_k), path="04_output_pathway_results/chicken_genetree_pathwayres_p0.05_q0.05.csv")

#Gene trees, all results, no cutoff
all_genes_nocutoff_k <- enrichKEGG(all_sig_gene,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=all_tested_gene,keyType="ncbi-geneid")

write_csv(summary(all_genes_nocutoff_k), path="04_output_pathway_results/chicken_genetree_pathwayres_nocutoffs.csv")

#Species trees, p<0.05, q<0.05
all_sp_k <- enrichKEGG(all_sig_sp,organism="gga",pvalueCutoff=0.05,pAdjustMethod="BH",qvalueCutoff=0.05,universe=all_tested_sp,keyType="ncbi-geneid")

write_csv(summary(all_sp_k), path="04_output_pathway_results/chicken_speciestree_pathwayres_p0.05_q0.05.csv")

#Species trees, all results, no cutoff

all_sp_nocutoff_k <- enrichKEGG(all_sig_sp,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=all_tested_sp,keyType="ncbi-geneid")

write_csv(summary(all_sp_nocutoff_k), path="04_output_pathway_results/chicken_speciestree_pathwayres_nocutoffs.csv")

#Plots *Will finalize when we decide what is going in the paper
#dotplot(all_genes_k)
#enrichMap(all_genes_k)
#cnetplot(all_genes_k,categorySize="pvalue",showCategory=10,fixed=TRUE)

#dotplot(all_sp_k)
#enrichMap(all_sp_k)
#cnetplot(all_sp_k,categorySize="pvalue",showCategory=10,fixed=TRUE)


#Extract proportion selected lineages (prop_sel.n) to visualize on pathways
prop_sel.n <- all_res_gene_ncbi %>%
  filter(!is.na(entrezgene)) %>%
  filter(prop_sel.n != "NaN") %>%
  pull(prop_sel.n)
names(prop_sel.n) <- all_res_gene_ncbi %>%
  filter(!is.na(entrezgene)) %>%
  filter(prop_sel.n != "NaN") %>%
  pull(entrezgene)


#pv.out.05164 <- pathview(gene.data=prop_sel.n,pathway.id="05164",species="gga",both.dirs=list(gene=F),mid="#e5f5f9",high="#2ca25f",limit=list(gene=max(prop_sel.n)),kegg.native=T,key.pos="topright",out.suffix="Influenza_A_propseln")
#pv.out.05164 <- pathview(gene.data=all_sig_gene,pathway.id="05164",species="gga",both.dirs=list(gene=F),mid="#e5f5f9",high="#2ca25f",kegg.native=T,key.pos="topright",out.suffix="Influenza_A_sig_genes")


prop_sel.n <- all_res_sp_ncbi %>%
  filter(!is.na(entrezgene)) %>%
  filter(prop_sel.n != "NaN") %>%
  pull(prop_sel.n)
names(prop_sel.n) <- all_res_sp_ncbi %>%
  filter(!is.na(entrezgene)) %>%
  filter(prop_sel.n != "NaN") %>%
  pull(entrezgene)


#pv.out.05164 <- pathview(gene.data=prop_sel.n,pathway.id="05164",species="gga",both.dirs=list(gene=F),mid="#e5f5f9",high="#2ca25f",limit=list(gene=max(prop_sel.n)),kegg.native=T,key.pos="topright",out.suffix="Influenza_A_propseln_sptree")
#pv.out.05164 <- pathview(gene.data=all_sig_sp,pathway.id="05164",species="gga",both.dirs=list(gene=F),mid="#e5f5f9",high="#2ca25f",kegg.native=T,key.pos="topright",out.suffix="Influenza_A_sig_genes_sptree")


#######################################################################################################################
#Pathway Enrichment Zebra Finch
#######################################################################################################################
load("02_output_annotated_data/all_res_zf_hs.Rdat")

#Testing with species tree
all_tested_sp_zf <- all_res_sp_zf_hs %>%
  filter(!is.na(FDRPval_m1m2) & !is.na(FDRPval_m2m2a) & !is.na(FDRPval_m7m8) & !is.na(FDRPval_m8m8a) & !is.na(FDRPval_busted)) %>%
  filter(!is.na(entrezgene_zf)) %>%
  distinct(hog,.keep_all=TRUE) %>%
  pull(entrezgene_zf)

all_sig_sp_zf <- all_res_sp_zf_hs %>%
  filter(!is.na(FDRPval_m1m2) & !is.na(FDRPval_m2m2a) & !is.na(FDRPval_m7m8) & !is.na(FDRPval_m8m8a) & !is.na(FDRPval_busted)) %>%
  filter(FDRPval_m1m2 < 0.05 & FDRPval_m2m2a < 0.05 & FDRPval_m7m8 < 0.05 & FDRPval_m8m8a < 0.05 & FDRPval_busted < 0.05) %>%
  filter(!is.na(entrezgene_zf)) %>%
  distinct(hog,.keep_all=TRUE) %>%
  pull(entrezgene_zf)

all_sp_zf_k <- enrichKEGG(all_sig_sp_zf,organism="tgu",pvalueCutoff=0.05,pAdjustMethod="BH",qvalueCutoff=0.05,universe=all_tested_sp_zf,keyType="ncbi-geneid")

write_csv(summary(all_sp_zf_k), path="04_output_pathway_results/zebrafinch_speciestree_pathwayres_p0.05_q0.05.csv")



#######################################################################################################################
#Pathway Enrichment Human
#######################################################################################################################

#Testing with species tree
all_tested_sp_hs <- all_res_sp_zf_hs %>%
  filter(!is.na(FDRPval_m1m2) & !is.na(FDRPval_m2m2a) & !is.na(FDRPval_m7m8) & !is.na(FDRPval_m8m8a) & !is.na(FDRPval_busted)) %>%
  filter(!is.na(entrezgene_hs)) %>%
  distinct(hog,.keep_all=TRUE) %>%
  pull(entrezgene_hs)

all_sig_sp_hs <- all_res_sp_zf_hs %>%
  filter(!is.na(FDRPval_m1m2) & !is.na(FDRPval_m2m2a) & !is.na(FDRPval_m7m8) & !is.na(FDRPval_m8m8a) & !is.na(FDRPval_busted)) %>%
  filter(FDRPval_m1m2 < 0.05 & FDRPval_m2m2a < 0.05 & FDRPval_m7m8 < 0.05 & FDRPval_m8m8a < 0.05 & FDRPval_busted < 0.05) %>%
  filter(!is.na(entrezgene_hs)) %>%
  distinct(hog,.keep_all=TRUE) %>%
  pull(entrezgene_hs)

all_sp_hs_k <- enrichKEGG(all_sig_sp_hs,organism="hsa",pvalueCutoff=0.05,pAdjustMethod="BH",qvalueCutoff=0.05,universe=all_tested_sp_hs,keyType="ncbi-geneid")


write_csv(summary(all_sp_hs_k), path="04_output_pathway_results/human_speciestree_pathwayres_p0.05_q0.05.csv")



#dotplot(all_sp_hs_k,showCategory=20)
#enrichMap(all_sp_hs_k)
#cnetplot(all_sp_hs_k,categorySize="pvalue",showCategory=20,fixed=TRUE)
