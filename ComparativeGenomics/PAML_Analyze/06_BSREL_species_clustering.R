setwd("~/Dropbox/BirdImmuneGeneEvolution")
library(tidyverse)
library(phytools)
library(RColorBrewer)
library(DOSE)
library(clusterProfiler)

#Load NCBI annotated dataset - output from script 02

load("02_output_annotated_data/all_res_ncbi.Rdat")
load("02_output_annotated_data/all_res_zf_hs.Rdat")



#######################################################################################################################
#Species clustering
#######################################################################################################################

#Using the abbreviations for all species, create a matrix of matrix of 1s and 0s for each hog, based on if that species' branch is significant (1), or not significant (1) - using the nominal set of branches. 
#Disregard all internal nodes.
#Note that in cases where species was missing from the hog dataset, we automatically assume it is not significant

species_info <- read_csv("06_input_cluster_by_species/species_list.csv")

#Remove outgroups
outgroups <- c("croPor","gavGan","anoCar","allMis","chrPic","cheMyd","allSin")
species_info <- species_info %>%
  filter(!(short_name %in% outgroups))

sp_abbr <- species_info$short_name

#Function to split string of sig branches/nodes, pull out and only return species names as a vector
format_branch_string <- function(branches){
  #split into relevant sections
  split_branches <- unlist(strsplit(branches,split = ":"))
  #remove nodes
  split_branches <- grep(pattern = "Node", x = split_branches,invert=TRUE,value=TRUE)
  #split by "_", only keep the first element (species abbreviation)
  split_branches <- purrr::map_chr(split_branches,function(x) unlist(strsplit(x,split="_"))[1])
  return(split_branches)
}

#Pull vector of nominally selected branch strings
nom_branches <- all_res_gene_ncbi %>%
  #dplyr::select(nom_branches) %>%
  mutate(nom_branches = as.character(nom_branches)) %>%
  pull(nom_branches)

#Presence/absence result matrix for all hogs
all_sel <- matrix(nrow=nrow(all_res_gene_ncbi), ncol=length(sp_abbr))

#Loop through all hogs and populate matrix
for (i in 1:length(nom_branches)){
  #Clean up names
  clean_branches <- format_branch_string(nom_branches[i])
  
  #Create presence/absence vector for all species that are in clean_branches
  all_sel[i,] <- sp_abbr %in% clean_branches
}
rownames(all_sel) <- all_res_gene_ncbi %>% pull(hog)
colnames(all_sel) <- sp_abbr

#PCA
all_sel_t <- t(all_sel)
pca_sp <- prcomp(all_sel_t)
summary(pca_sp)

loading <- pca_sp$rotation %>%
  as.data.frame %>%
  rownames_to_column("hog") %>%
  as.tibble

#Scree plot
pdf("06_output_cluster_by_species/pca_sp_scree_plot.pdf")
plot(pca_sp)
dev.off()


write_csv(data.frame(pca_sp$rotation),"06_output_cluster_by_species/pca_sp_loadings.csv")
write_csv(data.frame(pca_sp$x),"06_output_cluster_by_species/pca_sp_coordinates.csv")

pca_sp$x %>%
  as.data.frame %>%
  rownames_to_column() %>%
  as.tibble %>%
  ggplot(aes(PC1,PC2,label = rowname)) +
  geom_text()

sp_tree_bl <- read.tree("06_input_cluster_by_species/11700.tree1.nwk")
sp_tree <- read.tree("06_input_cluster_by_species/11700.final_spt.nwk")


pdf("06_output_cluster_by_species/PC_individual.pdf")
plotBranchbyTrait(sp_tree_bl,pca_sp$x[,"PC1"],mode="tips",palette = "heat.colors")
plotBranchbyTrait(sp_tree_bl,pca_sp$x[,"PC2"],mode="tips",palette = "heat.colors")
plotBranchbyTrait(sp_tree_bl,pca_sp$x[,"PC3"],mode="tips",palette = "heat.colors")
plotBranchbyTrait(sp_tree_bl,pca_sp$x[,"PC4"],mode="tips",palette = "heat.colors")
plotBranchbyTrait(sp_tree_bl,pca_sp$x[,"PC5"],mode="tips",palette = "heat.colors")
dev.off()

pdf("06_output_cluster_by_species/PCA_1.5_heatmap.pdf")
phylo.heatmap(tree=sp_tree,X=pca_sp$x[,1:6],colors = brewer.pal(n=9,name="PRGn"))
dev.off()
 

#Read in VIPs, BIPs and PIPs, see if any enrichment (or different PC scores) for different categories.
#Explore loadings
hog_anno <- read_csv("05_output_bird_mammal_comparison_results/hog_alt_annotation.csv") %>%
  mutate(hog = as.character(hog))

hog_geneids <- all_res_gene_zf_hs %>%
  dplyr::select(hog,entrezgene,entrezgene_zf,entrezgene_hs) %>%
  distinct(hog,.keep_all = TRUE)

loading_anno <- hog_geneids %>%
  left_join(hog_anno) %>%
  left_join(loading) %>%
  replace_na(list(vip = FALSE, bip = FALSE, pip = FALSE))

#Is there a difference in the distribution of PC scores for VIPs vs. non VIPs, etc?
loading_anno %>%
  ggplot(aes(PC1, fill=vip)) +
  geom_density(alpha=.5)
loading_anno %>%
  ggplot(aes(PC1, fill=bip)) +
  geom_density(alpha=.5)
loading_anno %>%
  ggplot(aes(PC1, fill=pip)) +
  geom_density(alpha=.5)


sp_coord <- pca_sp$x %>%
  as.data.frame %>%
  rownames_to_column("sp_abbr") %>%
  as.tibble

#Is there any enrichment in KEGG pathways, given PC loading scores?
PC1_loading_geneList <- loading_anno %>%
  filter(!is.na(entrezgene)) %>%
  pull(PC1) %>%
names(PC1_loading_geneList) <-loading_anno %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)
PC1_loading_geneList <- sort(PC1_loading_geneList,decreasing = TRUE)

pc1_kk <- gseKEGG(geneList = PC1_loading_geneList,organism = "gga", keyType = "ncbi-geneid", nPerm = 1000, pvalueCutoff = 1)
summary(pc1_kk)

PC2_loading_geneList <- loading_anno$PC2
names(PC2_loading_geneList) <- loading_anno$entrezgene
PC2_loading_geneList <- sort(PC2_loading_geneList,decreasing = TRUE)

pc2_kk <- gseKEGG(geneList = PC2_loading_geneList,organism = "gga", keyType = "ncbi-geneid", nPerm = 1000, pvalueCutoff = 1)
summary(pc2_kk)

PC3_loading_geneList <- loading_anno$PC3
names(PC3_loading_geneList) <- loading_anno$entrezgene
PC3_loading_geneList <- sort(PC3_loading_geneList,decreasing = TRUE)

pc3_kk <- gseKEGG(geneList = PC3_loading_geneList,organism = "gga", keyType = "ncbi-geneid", nPerm = 1000, pvalueCutoff = 1)
summary(pc3_kk)


#Enrichment in 5% tails of distribution?

all_tested <- loading_anno %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)

PC1_high_most_sig <- loading_anno %>%
  arrange(desc(PC1)) %>%
  head(n=0.025*nrow(loading_anno)) %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)
PC1_low_most_sig <- loading_anno %>%
  arrange(PC1) %>%
  head(n=0.025*nrow(loading_anno)) %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)
PC1_most_sig <- c(PC1_high_most_sig,PC1_low_most_sig)

PC1_k <- enrichKEGG(PC1_most_sig,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=all_tested,keyType="ncbi-geneid")
summary(PC1_k)


PC2_high_most_sig <- loading_anno %>%
  arrange(desc(PC2)) %>%
  head(n=0.025*nrow(loading_anno)) %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)
PC2_low_most_sig <- loading_anno %>%
  arrange(PC2) %>%
  head(n=0.025*nrow(loading_anno)) %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)
PC2_most_sig <- c(PC2_high_most_sig,PC2_low_most_sig)

PC2_k <- enrichKEGG(PC2_most_sig,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=all_tested,keyType="ncbi-geneid")
summary(PC2_k)


PC3_high_most_sig <- loading_anno %>%
  arrange(desc(PC3)) %>%
  head(n=0.025*nrow(loading_anno)) %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)
PC3_low_most_sig <- loading_anno %>%
  arrange(PC3) %>%
  head(n=0.025*nrow(loading_anno)) %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)
PC3_most_sig <- c(PC3_high_most_sig,PC3_low_most_sig)

PC3_k <- enrichKEGG(PC3_most_sig,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=all_tested,keyType="ncbi-geneid")
summary(PC3_k)

###############
#Combine with pathway enrichment results:
library(ggridges)

chicken_path_res <- read_csv("04_output_pathway_results/chicken_genetree_pathwayres_nocutoffs.csv")
chicken_path_res_ids <- chicken_path_res %>%
  filter(qvalue < 0.2) %>%
  separate_rows(geneID) %>%
  dplyr::select(Description,entrezgene = geneID)

sig_categories <- chicken_path_res_ids %>%
  distinct(Description) %>%
  pull(Description)

all_sig_gene <- all_res_gene_ncbi %>%
  filter(!is.na(FDRPval_m1m2) & !is.na(FDRPval_m2m2a) & !is.na(FDRPval_m7m8) & !is.na(FDRPval_m8m8a) & !is.na(FDRPval_busted)) %>%
  filter(!is.na(entrezgene)) %>%
  filter(FDRPval_m1m2 < 0.05 & FDRPval_m2m2a < 0.05 & FDRPval_m7m8 < 0.05 & FDRPval_m8m8a < 0.05 & FDRPval_busted < 0.05) %>%
  mutate(Description = "Sig genes") %>%
  dplyr::select(Description,entrezgene) %>%
  as.tibble

chicken_path_genes <- read_csv("04_output_pathway_results/chicken_genetree_allgenestested_pathwayres_nocutoffs.csv") %>%
  separate_rows(geneID) %>%
  dplyr::select(Description,entrezgene = geneID) 

all_geneIDs <- loading_anno %>%
  filter(!is.na(entrezgene)) %>%
  mutate(Description = "All genes") %>%
  dplyr::select(Description, entrezgene)

chicken_path_res_ids %>%
  bind_rows(all_geneIDs) %>%
  bind_rows(all_sig_gene) %>%
  left_join(loading_anno) %>%
  ggplot(aes(PC1,Description)) +
  geom_density_ridges()

chicken_path_genes %>%
  filter(Description %in% sig_categories) %>%
  left_join(loading_anno) %>%
  ggplot(aes(PC1,Description)) +
  geom_density_ridges()


loadings_classes <- chicken_path_res_ids %>%
  bind_rows(all_geneIDs) %>%
  bind_rows(all_sig_gene) %>%
  left_join(loading_anno)



all_genes_PC1 <- loadings_classes[loadings_classes$Description == "All genes", "PC1"] %>% pull
toll_PC1 <- loadings_classes[loadings_classes$Description == "Toll-like receptor signaling pathway", "PC1"] %>% pull
wilcox.test(all_genes_PC1,toll_PC1)

all_genes_PC1 <- loadings_classes[loadings_classes$Description == "All genes", "PC1"] %>% pull
hr_PC1 <- loadings_classes[loadings_classes$Description == "Homologous recombination", "PC1"] %>% pull
wilcox.test(all_genes_PC1,hr_PC1)


all_sig_genes_PC1 <- loadings_classes[loadings_classes$Description == "All genes" & loadings_classes$entrezgene %in% all_sig_gene, "PC2"] %>% pull
rand_1 <- sample(all_sig_genes_PC1, size=20, replace=F)
wilcox.test(all_genes_PC1,rand_1)
mean(all_sig_genes_PC1)
mean(rand_1)
mean(all_genes_PC1)
