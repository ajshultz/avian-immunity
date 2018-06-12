setwd("~/Dropbox/BirdImmuneGeneEvolution")
library(tidyverse)
library(phytools)
library(RColorBrewer)
library(DOSE)
library(clusterProfiler)
library(ape)
library(caper)
library(nlme)

#Load NCBI annotated dataset - output from script 02

load("02_output_annotated_data/all_res_ncbi.Rdat")
load("02_output_annotated_data/all_res_zf_hs.Rdat")



#######################################################################################################################
#Species clustering
#######################################################################################################################

#Using the abbreviations for all species, create a matrix of matrix of 1s and 0s for each hog, based on if that species' branch is significant (1), or not significant (1) - using the nominal set of branches. 
#Disregard all internal nodes.
#Note that in cases where species was missing from the hog dataset, we automatically assume it is not significant

species_info <- read_csv("06_input_cluster_by_species/sackton_et_al_species_list.csv")
lifespan_info <- read_csv("06_input_cluster_by_species/bird_lifespan_info.csv")

lifespan_info <- lifespan_info %>%
  dplyr::select(code,lifespan,captivity_or_wild)

#Remove outgroups
#outgroups <- c("croPor","gavGan","anoCar","allMis","chrPic","cheMyd","allSin")
#species_info <- species_info %>%
#  filter(!(short_name %in% outgroups))

sp_abbr <- species_info$code

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

#Get number of selected branches for each HOG
n_sig_sp <- all_sel %>%
  as.data.frame %>%
  rownames_to_column(var = "hog") %>%
  gather(sp_abbr,sel,-hog) %>%
  group_by(hog) %>%
  summarize(n_sel = sum(sel>0))


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

#####PCA with p-values as input
load("01_output_processed_data/bsrel_spval_res.RDat")


rep_na_mean <- function(vec){
  m <- mean(vec, na.rm = TRUE) 
  vec[is.na(vec)] <- m 
  return(vec) 
}

#How many missing hog results per species?
sp_missing_na <- bsrel_sp_pval_res_gene %>%
  filter(pval_type == "pval") %>%
  separate(hog_treenum,into=c("hog","treenum")) %>%
  dplyr::select(-pval_type,-treenum) %>%
  gather(sp,pval,-hog) %>%
  mutate(pval=ifelse(pval==0,1e-18,pval)) %>%
  mutate(pval=log10(pval)) %>%
  mutate(hog=paste0("hog_",hog)) %>%
  group_by(sp) %>%
  summarize(n_na = sum(is.na(pval)))

#Set all "0" values to 1e-18 (min is 1e-17), so can take log10
bsrel_sp_pval_res_gene_raw <- bsrel_sp_pval_res_gene %>%
  filter(pval_type == "pval") %>%
  separate(hog_treenum,into=c("hog","treenum")) %>%
  dplyr::select(-pval_type,-treenum) %>%
  gather(sp,pval,-hog) %>%
  mutate(pval=ifelse(pval==0,1e-18,pval)) %>%
  mutate(pval=log10(pval)) %>%
  mutate(hog=paste0("hog_",hog)) %>%
  spread(hog,pval) %>%
  mutate_at(vars(-sp),rep_na_mean)

#Which hogs are missing all data?
hogs_to_keep <- bsrel_sp_pval_res_gene %>%
  filter(pval_type == "pval") %>%
  separate(hog_treenum,into=c("hog","treenum")) %>%
  dplyr::select(-pval_type,-treenum) %>%
  gather(sp,pval,-hog) %>%
  mutate(pval=ifelse(pval==0,1e-18,pval)) %>%
  mutate(hog=paste0("hog_",hog)) %>%
  group_by(hog) %>%
  summarize(n_pvals = sum(!is.na(pval))) %>%
  filter(n_pvals > 0) %>%
  pull(hog)

#Prep format for PCA
bsrel_sp_pval_res_gene_forpca <- bsrel_sp_pval_res_gene_raw %>%
  dplyr::select(sp,hogs_to_keep) %>%
  as.data.frame
rownames(bsrel_sp_pval_res_gene_forpca) <- bsrel_sp_pval_res_gene_forpca$sp
bsrel_sp_pval_res_gene_forpca$sp <- NULL



pca_pvals <- prcomp(bsrel_sp_pval_res_gene_forpca)
summary(pca_pvals)

loading <- pca_pvals$rotation %>%
  as.data.frame %>%
  rownames_to_column("hog") %>%
  as.tibble

#Scree plot
pdf("06_output_cluster_by_species/pca_sp_pvals_scree_plot.pdf")
plot(pca_pvals)
dev.off()

write_csv(data.frame(pca_pvals$rotation),"06_output_cluster_by_species/pca_sp_pvals_loadings.csv")
write_csv(data.frame(pca_pvals$x),"06_output_cluster_by_species/pca_sp_pvals_coordinates.csv")

pca_pvals$x %>%
  as.data.frame %>%
  rownames_to_column() %>%
  as.tibble %>%
  ggplot(aes(PC1,PC2,label = rowname)) +
  geom_text()

pca_sp <- pca_pvals


#sp_tree_bl <- read.tree("06_input_cluster_by_species/11700.tree1.nwk")
sp_tree_bl <- read.tree("Alignments/10090.tree1.nwk")
sp_tree <- read.tree("06_input_cluster_by_species/11700.final_spt.nwk")


pdf("06_output_cluster_by_species/PC_individual.pdf")
plotBranchbyTrait(sp_tree_bl,pca_sp$x[,"PC1"],mode="tips",palette = "heat.colors")
plotBranchbyTrait(sp_tree_bl,pca_sp$x[,"PC2"],mode="tips",palette = "heat.colors")
plotBranchbyTrait(sp_tree_bl,pca_sp$x[,"PC3"],mode="tips",palette = "heat.colors")
plotBranchbyTrait(sp_tree_bl,pca_sp$x[,"PC4"],mode="tips",palette = "heat.colors")
plotBranchbyTrait(sp_tree_bl,pca_sp$x[,"PC5"],mode="tips",palette = "heat.colors")
dev.off()

pdf("06_output_cluster_by_species/PCA_1.5_heatmap.pdf")
phylo.heatmap(tree=sp_tree,X=pca_sp$x[,1:3],colors = brewer.pal(n=9,name="PRGn"))
dev.off()



 
#########Correlations with PC axes and life history traits


sp_coord <- pca_sp$x %>%
  as.data.frame %>%
  rownames_to_column("sp_abbr") %>%
  as.tibble

sp_coord_anno <- sp_coord %>%
  left_join(species_info,by=c("sp_abbr" = "code")) %>%
  mutate(heart_index = heart_mass_mean/body_mass_mean) %>%
  left_join(lifespan_info,by=c("sp_abbr" = "code"))

sp_coord_anno %>%
  ggplot(aes(PC1,log(body_mass_mean_hbabm))) +
  geom_point() +
  ylab("log(body mass)")
ggsave("06_output_cluster_by_species/PC1_log_BM_HABM.png")

sp_coord_anno %>%
  ggplot(aes(PC1,log(body_mass_mean))) +
  geom_point()

sp_coord_anno %>%
  ggplot(aes(lifespan,log(body_mass_mean_hbabm))) +
  geom_point()

sp_coord_anno %>%
  ggplot(aes(PC1,lifespan)) +
  geom_point()

sp_coord_anno %>%
  #filter(sp_abbr != "calAnn") %>%
  ggplot(aes(PC1,heart_index)) +
  geom_point()

sp_coord_anno %>%
  ggplot(aes(PC1,scaled_length_cr1)) +
  geom_point()

sp_coord_anno %>%
  ggplot(aes(PC1,assembly_size_gb)) +
  geom_point()

sp_coord_anno %>%
  ggplot(aes(PC1,agsd_genome_size_pg)) +
  geom_point()

sp_coord_anno %>%
  ggplot(aes(PC1,mass_specific_bmr)) +
  geom_point()

sp_coord_anno_df <- sp_coord_anno %>%
  as.data.frame


comp_dataset <- comparative.data(sp_tree_bl,sp_coord_anno_df, sp_abbr, vcv = T, vcv.dim = 3, na.omit=F)
hi_pc1 <- pgls(PC1 ~ heart_index, data = comp_dataset, lambda = "ML")
summary(hi_pc1)
hm_pc1 <- pgls(PC1 ~ log(heart_mass_median), data = comp_dataset, lambda = "ML")
summary(hm_pc1)
bm_pc1 <- pgls(PC1 ~ log(body_mass_median), data = comp_dataset, lambda = "ML")
summary(bm_pc1)
cr1_pc1 <- pgls(PC1 ~ scaled_length_cr1, data = comp_dataset, lambda = "ML")
summary(cr1_pc1)
bmr_pc1 <- pgls(PC1 ~ mass_specific_bmr, data = comp_dataset, lambda = "ML")
summary(bmr_pc1)


#Create reduced datasets with only species with body size/heart size data
red_dataset <- sp_coord_anno %>%
  #filter(sp_abbr != "calAnn") %>%
  filter(!is.na(heart_mass_mean)) %>%
  as.data.frame
rownames(red_dataset) <- red_dataset$sp_abbr

red_sp_tree_bl <- drop.tip(sp_tree_bl,tip = setdiff(sp_tree_bl$tip.label,red_dataset$sp_abbr))

#Run PGLS
comp_dataset <- comparative.data(red_sp_tree_bl,red_dataset, sp_abbr, vcv = T, vcv.dim = 3, na.omit=F)
hi_pc1 <- pgls(PC1 ~ heart_index, data = comp_dataset, lambda = "ML")
summary(hi_pc1)
hm_pc1 <- pgls(PC1 ~ log(heart_mass_mean), data = comp_dataset, lambda = "ML")
summary(hm_pc1)
bm_pc1 <- pgls(PC1 ~ log(body_mass_mean), data = comp_dataset, lambda = "ML")
summary(bm_pc1)



#Correlation with handbook of avian body masses
rownames(sp_coord_anno_df) <- sp_coord_anno_df$sp_abbr
bm_hbabm_pc1 <- gls(PC1 ~ log(body_mass_mean_hbabm),correlation=corBrownian(1,sp_tree_bl),data=sp_coord_anno_df)
summary(bm_hbabm_pc1)

#Correlation with lifespan
red_anno <- sp_coord_anno_df[!is.na(sp_coord_anno_df$lifespan),]
red_sp_tree_bl <- drop.tip(sp_tree_bl,setdiff(sp_tree_bl$tip.label,rownames(red_anno)))
bm_lifespan_pc1 <- gls(PC1 ~ lifespan,correlation=corBrownian(1,red_sp_tree_bl),data=red_anno)
summary(bm_lifespan_pc1)


#Phylogenetic regression of body size score and branch pval of each species

#Is there any enrichment in KEGG pathways, given PC loading scores?
#Pull vector of hogs that have at least 1 selected branch:
sel_hogs <- n_sig_sp %>% filter(n_sel > 1) %>% pull(hog)
names_sel_hogs <- paste0("hog_",sel_hogs)

hog_phylo_reg_res <- matrix(nrow=length(names_sel_hogs),ncol=3)
hog_phylo_reg_res[,1] <- sel_hogs

#Create a combo dataset with selected branches and species data
comb_data <- bsrel_sp_pval_res_gene_forpca %>%
  rownames_to_column(var = "sp_abbr") %>%
  as.tibble %>%
  left_join(sp_coord_anno)

bm_hog_res <- list()
#For each hog that has at least 1 selected branch, perform PGLS between mean body mass and selected branches. Record p-value, correlation
for (i in 1:length(names_sel_hogs)){
  hog <- names_sel_hogs[i]
  data <- comb_data %>%
    dplyr::select(sp_abbr,hog,body_mass_mean_hbabm) %>%
    as.data.frame
  rownames(data) <- data$sp_abbr
  
  hog_formula <- paste0(hog,"~log(body_mass_mean_hbabm)")
  
  try(bm_hog_res[[i]] <- gls(as.formula(hog_formula),correlation=corBrownian(1,sp_tree_bl),data=data))
  try(hog_phylo_reg_res[i,2] <- summary(bm_hog_res[[i]])$tTable[2,4])
  try(hog_phylo_reg_res[i,3] <- summary(bm_hog_res[[i]])$tTable[2,1])
}

colnames(hog_phylo_reg_res) <- c("hog","pvalue","coefficient")

hog_phylo_reg_res <- hog_phylo_reg_res %>%
  as.tibble %>%
  mutate(pvalue = as.double(pvalue), coefficient = as.double(coefficient))

save(bm_hog_res,hog_phylo_reg_res,file ="06_output_cluster_by_species/bodymass_sppvals_lm_res.Rdat")

hog_phylo_reg_res_anno <- hog_phylo_reg_res %>%
  left_join(all_res_gene_zf_hs) %>%
  dplyr::select(entrezgene,entrezgene_hs,pvalue,coefficient,FDRPval_busted,hog) %>%
  mutate(FDRPval_corr = p.adjust(pvalue,method="BH")) %>%
  arrange(FDRPval_corr)


comb_data %>%
  ggplot(aes(log(body_mass_mean_hbabm),sqrt(abs(hog_635)))) + geom_point()

comb_data %>%
  ggplot(aes(sqrt(abs(hog_635)))) + geom_histogram()






###Look at correlations with selected pvalues and body mass with Spearman's rank correlation, to account for non-normalcy (but doesn't account for phylogenetic structure)
#Is there any enrichment in KEGG pathways, given PC loading scores?
#Need to remake datasets to include missing data instead of using means, as we did for PCA
bsrel_sp_pval_res_gene_raw_na <- bsrel_sp_pval_res_gene %>%
  filter(pval_type == "pval") %>%
  separate(hog_treenum,into=c("hog","treenum")) %>%
  dplyr::select(-pval_type,-treenum) %>%
  gather(sp,pval,-hog) %>%
  mutate(pval=ifelse(pval==0,1e-18,pval)) %>%
  mutate(pval=log10(pval)) %>%
  mutate(hog=paste0("hog_",hog)) %>%
  spread(hog,pval)

#Which hogs are missing all data?
hogs_to_keep <- bsrel_sp_pval_res_gene %>%
  filter(pval_type == "pval") %>%
  separate(hog_treenum,into=c("hog","treenum")) %>%
  dplyr::select(-pval_type,-treenum) %>%
  gather(sp,pval,-hog) %>%
  mutate(pval=ifelse(pval==0,1e-18,pval)) %>%
  mutate(hog=paste0("hog_",hog)) %>%
  group_by(hog) %>%
  summarize(n_pvals = sum(!is.na(pval))) %>%
  filter(n_pvals > 0) %>%
  pull(hog)

#Prep format for PCA
bsrel_sp_pval_res_gene_forcorr <- bsrel_sp_pval_res_gene_raw_na %>%
  dplyr::select(sp,hogs_to_keep) %>%
  as.data.frame
rownames(bsrel_sp_pval_res_gene_forcorr) <- bsrel_sp_pval_res_gene_forcorr$sp
bsrel_sp_pval_res_gene_forcorr$sp <- NULL





#Pull vector of hogs that have at least 1 selected branch:
sel_hogs <- n_sig_sp %>% filter(n_sel > 1) %>% pull(hog)
names_sel_hogs <- paste0("hog_",sel_hogs)

hog_phylo_spear_res <- matrix(nrow=length(names_sel_hogs),ncol=4)
hog_phylo_spear_res[,1] <- sel_hogs

#Create a combo dataset with selected branches and species data
comb_data <- bsrel_sp_pval_res_gene_forcorr %>%
  rownames_to_column(var = "sp_abbr") %>%
  as.tibble %>%
  left_join(sp_coord_anno)

comb_data_df <- comb_data %>%
  #filter(!is.na(lifespan)) %>%
  as.data.frame

bm_hog_spear_res <- list()
#For each hog that has at least 1 selected branch, perform PGLS between mean body mass and selected branches. Record p-value, correlation
for (i in 1:length(names_sel_hogs)){
  hog <- names_sel_hogs[i]
  #data <- comb_data %>%
    #dplyr::select(sp_abbr,hog,body_mass_mean_hbabm) %>%
    #as.data.frame
  #rownames(data) <- data$sp_abbr
  
  x <- comb_data_df[,colnames(comb_data_df)==hog]
  y <- log(comb_data_df$body_mass_mean_hbabm)
  #y <- log(comb_data_df$lifespan)
  
  try(bm_hog_spear_res[[i]] <- cor.test(x,y,method="spearman"))
  try(hog_phylo_spear_res[i,2] <- bm_hog_spear_res[[i]]$p.value)
  try(hog_phylo_spear_res[i,3] <- bm_hog_spear_res[[i]]$estimate)
  try(hog_phylo_spear_res[i,4] <- bm_hog_spear_res[[i]]$statistic)
}

colnames(hog_phylo_spear_res) <- c("hog","pvalue","rho","statistic")

hog_phylo_spear_res <- hog_phylo_spear_res %>%
  as.tibble %>%
  mutate(pvalue = as.double(pvalue), rho = as.double(rho), statistic = as.double(statistic))

save(bm_hog_res,hog_phylo_spear_res,file ="06_output_cluster_by_species/bodymass_sppvals_spearman_res.Rdat")

hog_phylo_spear_res_anno <- hog_phylo_spear_res %>%
  left_join(all_res_gene_zf_hs) %>%
  dplyr::select(entrezgene,entrezgene_hs,pvalue,rho,FDRPval_busted,hog) %>%
  mutate(FDRPval_spear = p.adjust(pvalue,method="BH")) %>%
  left_join(n_sig_sp) %>%
  arrange(FDRPval_spear)


#hog_phylo_spear_res_anno <- hog_phylo_spear_res_anno %>%
#  left_join(n_sig_sp,by="hog")

hog_phylo_spear_res_anno %>%
  ggplot(aes(factor(n_sel),FDRPval_spear)) +
  geom_boxplot() +
  theme_bw()
ggsave("06_output_cluster_by_species/Spear_pval_by_nsp.pdf",width=10,height=4)






#Phylogenetic regression of body size score and presence or absense of each species
hog_phylo_reg_res <- matrix(nrow=length(names_sel_hogs),ncol=3)

rownames(all_sel) <- paste0("hog_",rownames(all_sel))
colnames(all_sel_t) <- paste0("hog_",colnames(all_sel_t))

#Is there any enrichment in KEGG pathways, given PC loading scores?
#Pull vector of hogs that have at least 1 selected branch:
sel_hogs <- n_sig_sp %>% filter(n_sel > 1) %>% pull(hog)
names_sel_hogs <- paste0("hog_",sel_hogs)

hog_phylo_reg_res[,1] <- sel_hogs
#Create a combo dataset with selected branches and species data
comb_data <- all_sel_t %>%
  as.data.frame %>%
  rownames_to_column(var = "sp_abbr") %>%
  as.tibble %>%
  left_join(sp_coord_anno)

bm_hog_res <- list()
#For each hog that has at least 1 selected branch, perform PGLS between mean body mass and selected branches. Record p-value, correlation
for (i in 1:length(names_sel_hogs)){
  hog <- names_sel_hogs[i]
  data <- comb_data %>%
    dplyr::select(sp_abbr,hog,body_mass_mean) %>%
    as.data.frame
  
  red_dataset <- comb_data %>%
    filter(!is.na(body_mass_mean)) %>%
    dplyr::select(sp_abbr,hog,body_mass_mean) %>%
    as.data.frame
  rownames(red_dataset) <- red_dataset$sp_abbr
  red_sp_tree_bl <- drop.tip(sp_tree_bl,tip = setdiff(sp_tree_bl$tip.label,red_dataset$sp_abbr))
  
  hog_formula <- paste0(hog,"~log(body_mass_mean)")
  
  try(bm_hog_res[[i]] <- gls(as.formula(hog_formula),correlation=corBrownian(1,red_sp_tree_bl),data=red_dataset))
  try(hog_phylo_reg_res[i,2] <- summary(bm_hog_res[[i]])$tTable[2,4])
  try(hog_phylo_reg_res[i,3] <- summary(bm_hog_res[[i]])$tTable[2,1])
}

colnames(hog_phylo_reg_res) <- c("hog","pvalue","coefficient")

hog_phylo_reg_res <- hog_phylo_reg_res %>%
  as.tibble %>%
  mutate(pvalue = as.double(pvalue), coefficient = as.double(coefficient))

save(bm_hog_res,hog_phylo_reg_res,file ="06_output_cluster_by_species/bodymass_selected_lm_res.Rdat")

hog_phylo_reg_res_anno <- hog_phylo_reg_res %>%
  left_join(all_res_gene_zf_hs) %>%
  dplyr::select(entrezgene,entrezgene_hs,pvalue,coefficient,FDRPval_busted,hog) %>%
  mutate(FDRPval_corr = p.adjust(pvalue,method="BH")) %>%
  arrange(FDRPval_corr)








#Read in VIPs, BIPs and PIPs, see if any enrichment (or different PC scores) for different categories.
#Explore loadings
hog_anno <- read_csv("05_output_bird_mammal_comparison_results/hog_alt_annotation.csv") %>%
  mutate(hog = as.character(hog)) %>%
  left_join(n_sig_sp) %>%
  mutate(hog = paste0("hog_",hog))

hog_geneids <- all_res_gene_zf_hs %>%
  dplyr::select(hog,entrezgene,entrezgene_zf,entrezgene_hs) %>%
  distinct(hog,.keep_all = TRUE) %>%
  mutate(hog = paste0("hog_",hog))

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






PC1_loading_geneList <- loading_anno %>%
  filter(!is.na(entrezgene)) %>%
  filter(hog %in% paste0("hog_",sel_hogs)) %>%
  #filter(n_sel > 5) %>%
  pull(PC1)
names(PC1_loading_geneList) <-loading_anno %>%
  filter(!is.na(entrezgene)) %>%
  filter(hog %in% paste0("hog_",sel_hogs)) %>%
  #filter(n_sel > 5) %>%
  pull(entrezgene)
PC1_loading_geneList <- sort(PC1_loading_geneList,decreasing = TRUE)

pc1_kk <- gseKEGG(geneList = PC1_loading_geneList,organism = "gga", keyType = "ncbi-geneid", nPerm = 1000, pvalueCutoff = 1,minGSSize = 10)
pc1_kk_summary <- summary(pc1_kk)



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
  head(n=0.1*nrow(loading_anno)) %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)
PC1_low_most_sig <- loading_anno %>%
  arrange(PC1) %>%
  head(n=0.1*nrow(loading_anno)) %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)
PC1_most_sig <- c(PC1_high_most_sig,PC1_low_most_sig)

PC1_high_k <- enrichKEGG(PC1_high_most_sig,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=0.2,universe=all_tested,keyType="ncbi-geneid")
summary(PC1_high_k)
PC1_low_k <- enrichKEGG(PC1_low_most_sig,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=0.2,universe=all_tested,keyType="ncbi-geneid")
summary(PC1_low_k)


PC2_high_most_sig <- loading_anno %>%
  arrange(desc(PC2)) %>%
  head(n=0.25*nrow(loading_anno)) %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)
PC2_low_most_sig <- loading_anno %>%
  arrange(PC2) %>%
  head(n=0.25*nrow(loading_anno)) %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)
PC2_most_sig <- c(PC2_high_most_sig,PC2_low_most_sig)

PC2_high_k <- enrichKEGG(PC2_high_most_sig,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=0.2,universe=all_tested,keyType="ncbi-geneid")
summary(PC2_high_k)
PC2_low_k <- enrichKEGG(PC2_low_most_sig,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=0.2,universe=all_tested,keyType="ncbi-geneid")
summary(PC2_low_k)


PC3_high_most_sig <- loading_anno %>%
  arrange(desc(PC3)) %>%
  head(n=0.25*nrow(loading_anno)) %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)
PC3_low_most_sig <- loading_anno %>%
  arrange(PC3) %>%
  head(n=0.25*nrow(loading_anno)) %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)
PC3_most_sig <- c(PC3_high_most_sig,PC3_low_most_sig)

PC3_high_k <- enrichKEGG(PC3_high_most_sig,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=0.2,universe=all_tested,keyType="ncbi-geneid")
summary(PC3_high_k)
PC3_low_k <- enrichKEGG(PC3_low_most_sig,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=0.2,universe=all_tested,keyType="ncbi-geneid")
summary(PC3_low_k)





#GSE with correlation

Corr_loading_geneList <- hog_phylo_reg_res_anno %>%
  filter(!is.na(entrezgene)) %>%
  pull(coefficient)
names(Corr_loading_geneList) <-hog_phylo_reg_res_anno %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)
Corr_loading_geneList <- sort(Corr_loading_geneList,decreasing = TRUE)

corr_kk <- gseKEGG(geneList = Corr_loading_geneList,organism = "gga", keyType = "ncbi-geneid", nPerm = 1000, pvalueCutoff = 1,minGSSize = 10)
test <- summary(corr_kk)

#Significant pvalues?
corr_sig <- hog_phylo_reg_res_anno %>%
  filter(FDRPval_corr < 0.05) %>%
  pull(entrezgene)
all_tested <- hog_phylo_reg_res_anno %>%
  pull(entrezgene)
Corr_k <- enrichKEGG(corr_sig,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=0.2,universe=all_tested,keyType="ncbi-geneid")
summary(Corr_k)



#GSE with spearman correlation

Corr_loading_geneList <- hog_phylo_spear_res_anno %>%
  filter(!is.na(entrezgene)) %>%
  filter(n_sel > 5) %>%
  pull(rho)
names(Corr_loading_geneList) <-hog_phylo_spear_res_anno %>%
  filter(!is.na(entrezgene)) %>%
  filter(n_sel > 5) %>%
  pull(entrezgene)
Corr_loading_geneList <- sort(Corr_loading_geneList,decreasing = TRUE)

corr_kk <- gseKEGG(geneList = Corr_loading_geneList,organism = "gga", keyType = "ncbi-geneid", nPerm = 1000, pvalueCutoff = 1,minGSSize = 10,by="DOSE")
corr_kk_summary <- summary(corr_kk)

#Significant pvalues?
corr_sig <- hog_phylo_spear_res_anno %>%
  filter(FDRPval_spear < 0.2) %>%
  filter(!is.na(entrezgene)) %>%
  filter(n_sel > 5) %>%
  pull(entrezgene)
all_tested <- hog_phylo_spear_res_anno %>%
  filter(n_sel > 5) %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)
Corr_k <- enrichKEGG(corr_sig,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=all_tested,keyType="ncbi-geneid")
spear_enrich <- summary(Corr_k)
spear_enrich[1:5,]

#Plot cellular senescence pathway
library(pathview)
rho <- hog_phylo_spear_res_anno %>% pull(rho)
names(rho) <- hog_phylo_spear_res_anno %>% pull(entrezgene)
pv.out.04218 <- pathview(gene.data=rho,pathway.id="04218",species="gga",both.dirs=list(gene=F),mid="#e5f5f9",high="#2ca25f",limit=list(gene=max(rho)),kegg.native=T,key.pos="topright",out.suffix="Cellular_senescence_rho")



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
  bind_rows(all_geneIDs) %>%
  filter(Description %in% sig_categories | Description == "All genes") %>%
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



