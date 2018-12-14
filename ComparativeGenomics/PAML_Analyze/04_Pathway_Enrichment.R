library(tidyverse)
library(DOSE)
library(clusterProfiler)
library(pathview)
library(forcats)
library(cowplot)

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


#######################################################################################################################
#Pathway Enrichment
#######################################################################################################################

#Chicken pathway enrichemnt tests:
#Gene trees, q<0.1
all_genes_k <- enrichKEGG(all_sig_gene,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=0.1,universe=all_tested_gene,keyType="ncbi-geneid")

write_csv(as.data.frame(all_genes_k), path="04_output_pathway_results/chicken_genetree_pathwayres_p1_q0.1.csv")

as.data.frame(all_genes_k)

#Gene trees, all results, no cutoff
all_genes_nocutoff_k <- enrichKEGG(all_sig_gene,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=all_tested_gene,keyType="ncbi-geneid")

write_csv(as.data.frame(all_genes_nocutoff_k), path="04_output_pathway_results/chicken_genetree_pathwayres_nocutoffs.csv")

#Which genes are present in each pathway for this dataset?
all_genes_tested_nocutoff_k <- enrichKEGG(all_tested_gene,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=all_tested_gene,keyType="ncbi-geneid")

write_csv(as.data.frame(all_genes_tested_nocutoff_k), path="04_output_pathway_results/chicken_genetree_allgenestested_pathwayres_nocutoffs.csv")

#Species trees, q<0.1
all_sp_k <- enrichKEGG(all_sig_sp,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=0.1,universe=all_tested_sp,keyType="ncbi-geneid")

write_csv(as.data.frame(all_sp_k), path="04_output_pathway_results/chicken_speciestree_pathwayres_p1_q0.1.csv")

#Species trees, all results, no cutoff

all_sp_nocutoff_k <- enrichKEGG(all_sig_sp,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=all_tested_sp,keyType="ncbi-geneid")

write_csv(as.data.frame(all_sp_nocutoff_k), path="04_output_pathway_results/chicken_speciestree_pathwayres_nocutoffs.csv")

#Plots *Will finalize when we decide what is going in the paper
#dotplot(all_genes_k,showCategory=30)

## count the gene number
all_genes_dotplot <- all_genes_k %>%
  as.tibble %>%
  separate(GeneRatio,into=c("sig_genes_pathway","sig_genes"),remove = F) %>%
  separate(BgRatio, into=c("bg_genes_pathway","bg_genes"),remove = F) %>%
  mutate(GeneRatio = as.numeric(sig_genes_pathway)/as.numeric(sig_genes)) %>%
  mutate(enrichment = (as.numeric(sig_genes_pathway)/as.numeric(sig_genes))/(as.numeric(bg_genes_pathway)/as.numeric(bg_genes))) %>%
  ggplot(aes(x = enrichment, y = fct_reorder(Description, enrichment))) + 
  geom_point(aes(size = enrichment, color = qvalue)) +
  theme_bw(base_size = 12) +
  scale_color_gradient2(high="#332288", mid = "#DDCC77") +
  guides(colour = guide_colorbar(reverse=T)) +
  ylab(NULL) +
  xlab("enrichment") +
  theme(axis.text=element_text(color="black"))

all_genes_dotplot
ggsave("04_output_pathway_results/chicken_genetree_qval0.1_dotplot.pdf",width=6,height=5)


#Produce cnetplot to show how genes from different pathways overlap
all_genes_cnet <- cnetplot(all_genes_k,showCategory = 30)

#Need to report how many sig categories to show, so that we can turn off gene names
n_sig_categories <- 18

#Change gene names to blank spaces
all_genes_cnet$data$name <- c(as.character(all_genes_cnet$data$name[1:n_sig_categories]),rep("",nrow(all_genes_cnet$data)-n_sig_categories))

#Change colors to desired hex codes
node_col <- "#332288"
gene_col <- "#44AA99"
all_genes_cnet$data$color <- c(rep(node_col,n_sig_categories),rep(gene_col,nrow(all_genes_cnet$data)-n_sig_categories))

#Reduce size of genes a bit
all_genes_cnet$data$size <- c(rep(-1,n_sig_categories),rep(-4,nrow(all_genes_cnet$data)-n_sig_categories))
all_genes_cnet$data$size <- c(all_genes_cnet$data$size[1:n_sig_categories],rep(-4,nrow(all_genes_cnet$data)-n_sig_categories))

#Reduce text size
all_genes_cnet$theme$text$size=10

#Save
all_genes_cnet
ggsave("04_output_pathway_results/chicken_genetree_qval0.1_cnetplot.pdf",width=8,height=8)

#Plot both together and save
plot_grid(all_genes_dotplot,all_genes_cnet,ncol=1,labels=c("A","B"),scale = 0.9)
ggsave("04_output_pathway_results/chicken_genetree_dotplot_cnet_figure.pdf",width=6,height=10)

#dotplot(all_sp_k)
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

#Create directory for pathway figs if it doesn't exit
if (!("04_output_pathway_results/pathway_figs" %in% list.dirs())){
  dir.create("04_output_pathway_results/pathway_figs")
}

#Have to re-set working direcotry so figures are produced in correct location
setwd("04_output_pathway_results/pathway_figs")
pv.out.05164 <- pathview(gene.data=prop_sel.n,pathway.id=all_genes_k$ID,species="gga",both.dirs=list(gene=F),mid="#e5f5f9",high="#2ca25f",limit=list(gene=max(prop_sel.n)),kegg.native=T,key.pos="topright",out.suffix="pathway_genetrees_propseln")
pv.out.05164 <- pathview(gene.data=all_sig_gene,pathway.id=all_genes_k$ID,species="gga",both.dirs=list(gene=F),mid="#e5f5f9",high="#2ca25f",kegg.native=T,key.pos="topright",out.suffix="pathway_genetrees_sig_genes")


prop_sel.n <- all_res_sp_ncbi %>%
  filter(!is.na(entrezgene)) %>%
  filter(prop_sel.n != "NaN") %>%
  pull(prop_sel.n)
names(prop_sel.n) <- all_res_sp_ncbi %>%
  filter(!is.na(entrezgene)) %>%
  filter(prop_sel.n != "NaN") %>%
  pull(entrezgene)


pv.out.05164 <- pathview(gene.data=prop_sel.n,pathway.id=all_sp_k$ID,species="gga",both.dirs=list(gene=F),mid="#e5f5f9",high="#2ca25f",limit=list(gene=max(prop_sel.n)),kegg.native=T,key.pos="topright",out.suffix="pathway_sptrees_propseln")
pv.out.05164 <- pathview(gene.data=all_sig_sp,pathway.id=all_sp_k$ID,species="gga",both.dirs=list(gene=F),mid="#e5f5f9",high="#2ca25f",kegg.native=T,key.pos="topright",out.suffix="pathway_sptrees_sig_genes")

setwd("../../")

#######################################################################################################################
#Pathway Enrichment Zebra Finch
#######################################################################################################################
load("02_output_annotated_data/all_res_zf_hs.Rdat")

#Testing with gene tree
all_tested_gene_zf <- all_res_gene_zf_hs %>%
  filter(!is.na(FDRPval_m1m2) & !is.na(FDRPval_m2m2a) & !is.na(FDRPval_m7m8) & !is.na(FDRPval_m8m8a) & !is.na(FDRPval_busted)) %>%
  filter(!is.na(entrezgene_zf)) %>%
  distinct(hog,.keep_all=TRUE) %>%
  pull(entrezgene_zf)

all_sig_gene_zf <- all_res_gene_zf_hs %>%
  filter(!is.na(FDRPval_m1m2) & !is.na(FDRPval_m2m2a) & !is.na(FDRPval_m7m8) & !is.na(FDRPval_m8m8a) & !is.na(FDRPval_busted)) %>%
  filter(FDRPval_m1m2 < 0.05 & FDRPval_m2m2a < 0.05 & FDRPval_m7m8 < 0.05 & FDRPval_m8m8a < 0.05 & FDRPval_busted < 0.05) %>%
  filter(!is.na(entrezgene_zf)) %>%
  distinct(hog,.keep_all=TRUE) %>%
  pull(entrezgene_zf)

all_gene_zf_k <- enrichKEGG(all_sig_gene_zf,organism="tgu",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=all_tested_gene_zf,keyType="ncbi-geneid")

write_csv(as.data.frame(all_gene_zf_k), path="04_output_pathway_results/zebrafinch_genetree_pathwayres_nocutoff.csv")

all_gene_zf_k_q0.1 <- enrichKEGG(all_sig_gene_zf,organism="tgu",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=0.1,universe=all_tested_gene_zf,keyType="ncbi-geneid")


#Produce cnetplot to show how genes from different pathways overlap
all_genes_cnet <- cnetplot(all_gene_zf_k_q0.1,showCategory = 30)

#Need to report how many sig categories to show, so that we can turn off gene names
n_sig_categories <- nrow(as.data.frame(all_gene_zf_k_q0.1))

#Change gene names to blank spaces
all_genes_cnet$data$name <- c(as.character(all_genes_cnet$data$name[1:n_sig_categories]),rep("",nrow(all_genes_cnet$data)-n_sig_categories))

#Change colors to desired hex codes
node_col <- "#332288"
gene_col <- "#44AA99"
all_genes_cnet$data$color <- c(rep(node_col,n_sig_categories),rep(gene_col,nrow(all_genes_cnet$data)-n_sig_categories))

#Reduce size of genes a bit
all_genes_cnet$data$size <- c(rep(-1,n_sig_categories),rep(-4,nrow(all_genes_cnet$data)-n_sig_categories))
all_genes_cnet$data$size <- c(all_genes_cnet$data$size[1:n_sig_categories],rep(-4,nrow(all_genes_cnet$data)-n_sig_categories))

#Reduce text size
all_genes_cnet$theme$text$size=10

#Save
all_genes_cnet
ggsave("04_output_pathway_results/zebrafinch_genetree_qval0.1_cnetplot.pdf",width=8,height=8)


#######################################################################################################################
#Pathway Enrichment Human
#######################################################################################################################

#Testing with species tree
all_tested_gene_hs <- all_res_gene_zf_hs %>%
  filter(!is.na(FDRPval_m1m2) & !is.na(FDRPval_m2m2a) & !is.na(FDRPval_m7m8) & !is.na(FDRPval_m8m8a) & !is.na(FDRPval_busted)) %>%
  filter(!is.na(entrezgene_hs)) %>%
  distinct(hog,.keep_all=TRUE) %>%
  pull(entrezgene_hs)

all_sig_gene_hs <- all_res_gene_zf_hs %>%
  filter(!is.na(FDRPval_m1m2) & !is.na(FDRPval_m2m2a) & !is.na(FDRPval_m7m8) & !is.na(FDRPval_m8m8a) & !is.na(FDRPval_busted)) %>%
  filter(FDRPval_m1m2 < 0.05 & FDRPval_m2m2a < 0.05 & FDRPval_m7m8 < 0.05 & FDRPval_m8m8a < 0.05 & FDRPval_busted < 0.05) %>%
  filter(!is.na(entrezgene_hs)) %>%
  distinct(hog,.keep_all=TRUE) %>%
  pull(entrezgene_hs)

all_gene_hs_k <- enrichKEGG(all_sig_gene_hs,organism="hsa",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=all_tested_gene_hs,keyType="ncbi-geneid")

write_csv(as.data.frame(all_gene_hs_k), path="04_output_pathway_results/human_speciestree_pathwayres_nocutoff.csv")

all_gene_hs_k_q0.1 <- enrichKEGG(all_sig_gene_hs,organism="hsa",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=0.1,universe=all_tested_gene_hs,keyType="ncbi-geneid")

#Produce cnetplot to show how genes from different pathways overlap
all_genes_cnet <- cnetplot(all_gene_hs_k_q0.1,showCategory = 40)

#Need to report how many sig categories to show, so that we can turn off gene names
n_sig_categories <- nrow(as.data.frame(all_gene_hs_k_q0.1))

#Change gene names to blank spaces
all_genes_cnet$data$name <- c(as.character(all_genes_cnet$data$name[1:n_sig_categories]),rep("",nrow(all_genes_cnet$data)-n_sig_categories))

#Change colors to desired hex codes
node_col <- "#332288"
gene_col <- "#44AA99"
all_genes_cnet$data$color <- c(rep(node_col,n_sig_categories),rep(gene_col,nrow(all_genes_cnet$data)-n_sig_categories))

#Reduce size of genes a bit
all_genes_cnet$data$size <- c(rep(-1,n_sig_categories),rep(-4,nrow(all_genes_cnet$data)-n_sig_categories))
all_genes_cnet$data$size <- c(all_genes_cnet$data$size[1:n_sig_categories],rep(-4,nrow(all_genes_cnet$data)-n_sig_categories))

#Reduce text size
all_genes_cnet$theme$text$size=10

#Save
all_genes_cnet
ggsave("04_output_pathway_results/human_genetree_qval0.1_cnetplot.pdf",width=12,height=12)



#Is there any difference in gene lengths by pathway?

#Use all tested genes as enrichment set, to use gene output from clusterProfiler to grab all genes belonging in a given pathway
all_tested_gene <- all_res_gene_ncbi %>%
  filter(!is.na(FDRPval_m1m2) & !is.na(FDRPval_m2m2a) & !is.na(FDRPval_m7m8) & !is.na(FDRPval_m8m8a) & !is.na(FDRPval_busted)) %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)

all_tested_sp <- all_res_sp_ncbi %>%
  filter(!is.na(FDRPval_m1m2) & !is.na(FDRPval_m2m2a) & !is.na(FDRPval_m7m8) & !is.na(FDRPval_m8m8a) & !is.na(FDRPval_busted)) %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)

#Chicken pathway enrichemnt tests:
all_genes_testing <- enrichKEGG(all_tested_gene,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=all_tested_gene,keyType="ncbi-geneid")

all_genes_pathways_lengths <- all_genes_testing %>% as.data.frame %>%
  as.tibble %>%
  separate_rows(geneID) %>%
  left_join(all_res_gene_ncbi,by=c("geneID" = "entrezgene")) %>%
  dplyr::select(ID,Description,geneID,length)
  
sig_pathways <- all_genes_k$Description

all_genes_pathways_lengths_median <- all_genes_pathways_lengths %>%
  group_by(Description) %>%
  summarize(median_length = median(length),sd_length = sd(length)) %>%
  arrange(desc(median_length)) %>%
  mutate(sig_pathway = if_else(Description %in% sig_pathways,TRUE,FALSE)) %>%
  print(n=150)

write_csv(all_genes_pathways_lengths_median,"04_output_pathway_results/chicken_genetree_pathway_lengths.csv")

#Test for differences in median length between pathways identified as being significantly enriched in the "immune and signaling" or "recombination and DNA repair" KEGG categories vs. all other categories.
immune <- c("Cytokine-cytokine receptor interaction","Influenza A","Necroptosis","ECM-receptor interaction","Herpes simplex infection", "Toll-like receptor signaling pathway","Apoptosis","Phagosome","Cell adhesion molecules (CAMs)","RIG-I-like receptor signaling pathway")
recomb <- c("Fanconi anemia pathway","Homologous recombination","Base excision repair","Non-homologous end-joining","Mismatch repair")

all_genes_pathways_lengths_median <- all_genes_pathways_lengths_median %>%
  mutate(immune = if_else(Description %in% immune,TRUE,FALSE),
         recomb = if_else(Description %in% recomb,TRUE,FALSE))


sink("04_output_pathway_results/chicken_genetree_mann_whitney_for_alignment_length.txt")
#Immune pathways
with(all_genes_pathways_lengths_median,wilcox.test(median_length~immune))

#Recombination pathways
with(all_genes_pathways_lengths_median,wilcox.test(median_length~recomb))
sink()


