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

species_info <- read_csv("06_input_cluster_by_species/sackton_et_al_species_list.csv")

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


#####PCA with p-values as input
load("01_output_processed_data/bsrel_spval_res.RDat")

#Replace NA values with mean for hog
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
  mutate(hog=paste0("hog_",hog))

#Extract omega values (for second rate class)
bsrel_sp_omega_res_gene_raw <- bsrel_sp_pval_res_gene %>%
  filter(pval_type == "omega") %>%
  separate(hog_treenum,into=c("hog","treenum")) %>%
  dplyr::select(-pval_type,-treenum) %>%
  gather(sp,omega,-hog) %>%
  mutate(hog=paste0("hog_",hog))

#Extract nubmer rate classes values (for second rate class)
bsrel_sp_rate_classes_res_gene_raw <- bsrel_sp_pval_res_gene %>%
  filter(pval_type == "rate_classes") %>%
  separate(hog_treenum,into=c("hog","treenum")) %>%
  dplyr::select(-pval_type,-treenum) %>%
  gather(sp,rate_classes,-hog) %>%
  mutate(hog=paste0("hog_",hog))

#Pull out the weights (proportion of sites under selection)
#Set weights of genes not identified as being under selection to 0 (pval > 0.05) and weights of genes identified as having one rate class to 0, do the same with omega. Combene all parameters into one tibble
#Set cap of omega value to 100
bsrel_sp_all_params_res_gene_raw <- bsrel_sp_pval_res_gene %>%
  filter(pval_type == "weight") %>%
  separate(hog_treenum,into=c("hog","treenum")) %>%
  dplyr::select(-pval_type,-treenum) %>%
  gather(sp,weight,-hog) %>%
  mutate(hog=paste0("hog_",hog)) %>%
  left_join(bsrel_sp_pval_res_gene_raw) %>%
  left_join(bsrel_sp_omega_res_gene_raw) %>%
  left_join(bsrel_sp_rate_classes_res_gene_raw) %>%
  mutate(weight_sig = if_else(rate_classes == 1,0,weight)) %>%
  mutate(weight_sig = if_else(pval > 0.05,0,weight_sig)) %>%
  mutate(omega_sig = if_else(rate_classes == 1,0,omega)) %>%
  mutate(omega_sig = if_else(pval > 0.05,0,omega_sig)) %>%
  mutate(omega_sig = if_else(omega_sig > 10000,10000,omega_sig)) %>%
  mutate(omega_weight_sig = omega_sig * weight_sig)

bsrel_sp_all_params_res_gene_raw %>%
  filter(omega_sig > 0) %>%
  ggplot(aes(log(omega_sig))) +
  geom_histogram()

bsrel_sp_all_params_res_gene_raw %>%
  filter(weight_sig > 0) %>%
  ggplot(aes(weight_sig)) +
  geom_histogram()

bsrel_sp_all_params_res_gene_raw %>%
  filter(omega_sig > 0, weight_sig > 0) %>%
  ggplot(aes(log(omega_weight_sig),log(pval))) +
  geom_point(alpha=0.01)

bsrel_sp_all_params_res_gene_raw_trunc <- bsrel_sp_all_params_res_gene_raw %>%
  filter(omega_sig > 0, weight_sig > 0,pval!=1e-18)
summary(lm(log(omega_weight_sig)~log(pval),data=bsrel_sp_all_params_res_gene_raw_trunc))
  
bsrel_sp_weight_res_gene_raw <- bsrel_sp_all_params_res_gene_raw %>%
  dplyr::select(sp,hog,weight_sig) %>%
  spread(hog,weight_sig) %>%
  mutate_at(vars(-sp),rep_na_mean)

bsrel_sp_omega_res_gene_raw <- bsrel_sp_all_params_res_gene_raw %>%
  dplyr::select(sp,hog,omega_sig) %>%
  #mutate(omega_sig=ifelse(omega_sig==0,1e-18,omega_sig)) %>%
  spread(hog,omega_sig) %>%
  mutate_at(vars(-sp),rep_na_mean)


log_use_0s_omega <- function(column){
  column[column == 0] <- 1
  return(log(column))
}

#For omega, replace all 0 values with 1 (neutral value), then log transform
bsrel_sp_omega_res_gene_trans <- bsrel_sp_omega_res_gene_raw %>%
  mutate_at(vars(-sp),log_use_0s_omega)

bsrel_sp_pval_res_gene_raw <- bsrel_sp_all_params_res_gene_raw %>%
  dplyr::select(sp,hog,pval) %>%
  mutate(pval=log10(pval)) %>%
  spread(hog,pval)  %>%
  mutate_at(vars(-sp),rep_na_mean)


save(bsrel_sp_all_params_res_gene_raw, bsrel_sp_omega_res_gene_trans, bsrel_sp_weight_res_gene_raw, bsrel_sp_pval_res_gene_raw, file = "06_output_cluster_by_species/parameter_datasets_transformed.Rdat")
#################Playing around with plotting omega and weight
bsrel_sp_all_params_res_gene_raw %>%
  arrange(desc(omega_sig))

bsrel_sp_all_params_res_gene_raw %>%
  filter(omega_sig == 1e26) %>%
  mutate(omega_sig = ifelse(omega_sig == 0, 1, omega_sig)) %>%
  ggplot(aes(log(pval),log(omega_sig))) +
  geom_bin2d()

bsrel_sp_all_params_res_gene_raw %>%
  filter(omega_sig != 1e26) %>%
  filter(omega_sig != 0) %>%
  ggplot(aes(log(pval),log(omega_sig))) +
  geom_point(alpha=0.01)

bsrel_sp_all_params_res_gene_raw %>%
  filter(weight_sig != 0) %>%
  ggplot(aes(log(pval),weight_sig)) +
  geom_point(alpha=0.01)

bsrel_sp_all_params_res_gene_raw %>%
  filter(omega < 9000) %>%
  mutate(weight_omega_sig = weight_sig * omega_sig) %>%
  filter(weight_omega_sig != 0) %>%
  ggplot(aes(log(pval),weight_omega_sig)) +
  geom_point(alpha=0.01)
#########################

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

bsrel_sp_weight_res_gene_forpca <- bsrel_sp_weight_res_gene_raw %>%
  dplyr::select(sp,hogs_to_keep) %>%
  as.data.frame
rownames(bsrel_sp_weight_res_gene_forpca) <- bsrel_sp_weight_res_gene_forpca$sp
bsrel_sp_weight_res_gene_forpca$sp <- NULL

bsrel_sp_omega_res_gene_forpca <- bsrel_sp_omega_res_gene_trans %>%
  dplyr::select(sp,hogs_to_keep) %>%
  as.data.frame
rownames(bsrel_sp_omega_res_gene_forpca) <- bsrel_sp_omega_res_gene_forpca$sp
bsrel_sp_omega_res_gene_forpca$sp <- NULL


#Run PCA
pca_pvals <- prcomp(bsrel_sp_pval_res_gene_forpca)
pca_omega <- prcomp(bsrel_sp_omega_res_gene_forpca)
pca_weights <- prcomp(bsrel_sp_weight_res_gene_forpca)

#Extract loadings
loading <- pca_pvals$rotation %>%
  as.data.frame %>%
  rownames_to_column("hog") %>%
  as.tibble

loading_weights <- pca_weights$rotation %>%
  as.data.frame %>%
  rownames_to_column("hog") %>%
  as.tibble

loading_omega <- pca_omega$rotation %>%
  as.data.frame %>%
  rownames_to_column("hog") %>%
  as.tibble

#Scree plot
pdf("06_output_cluster_by_species/pca_gene_pvals_scree_plot.pdf")
plot(pca_pvals)
dev.off()

pdf("06_output_cluster_by_species/pca_gene_omega_scree_plot.pdf")
plot(pca_omega)
dev.off()

pdf("06_output_cluster_by_species/pca_gene_weights_scree_plot.pdf")
plot(pca_weights)
dev.off()

#Variance
summary(pca_pvals)
summary(pca_omega)
summary(pca_weights)

#Save pval results
write_csv(data.frame(pca_pvals$rotation),"06_output_cluster_by_species/pca_gene_pvals_loadings.csv")
write_csv(data.frame(pca_pvals$x),"06_output_cluster_by_species/pca_gene_pvals_coordinates.csv")

save(pca_pvals,bsrel_sp_pval_res_gene_raw, file="06_output_cluster_by_species/pca_gene_pvals_all_res.rDat")

#Save omega results
write_csv(data.frame(pca_omega$rotation),"06_output_cluster_by_species/pca_gene_omega_loadings.csv")
write_csv(data.frame(pca_omega$x),"06_output_cluster_by_species/pca_gene_omega_coordinates.csv")

save(pca_omega,bsrel_sp_omega_res_gene_raw,bsrel_sp_omega_res_gene_trans, file="06_output_cluster_by_species/pca_gene_omega_all_res.rDat")

#Save weights results
write_csv(data.frame(pca_weights$rotation),"06_output_cluster_by_species/pca_gene_weights_loadings.csv")
write_csv(data.frame(pca_weights$x),"06_output_cluster_by_species/pca_gene_weights_coordinates.csv")

save(pca_weights,bsrel_gene_weight_res_gene_raw, file="06_output_cluster_by_species/pca_gene_weights_all_res.rDat")







########Do all anlayses using pvalues
pca_sp <- pca_pvals

#Read in tree, chose random hog to give branch lengths. Species tree has same topology for all hogs. Note that results are robust to hog selected.
sp_tree_bl <- read.tree("06_input_cluster_by_species/11700.tree1.nwk")
#sp_tree_bl <- read.tree("Alignments/10090.tree1.nwk")
sp_tree <- read.tree("06_input_cluster_by_species/11700.final_spt.nwk")

#Making custom version of plotBranchbyTrait from phytools to include plotting without branch lengths and custom color palette.
plotBranchbyTrait<-function(tree,x,mode=c("edges","tips","nodes"),palette="rainbow",legend=TRUE,xlims=NULL,...){
  mode<-mode[1]
  if(mode=="tips"){
    x<-c(x[tree$tip.label],fastAnc(tree,x))
    names(x)[1:length(tree$tip.label)]<-1:length(tree$tip.label)
    XX<-matrix(x[tree$edge],nrow(tree$edge),2)
    x<-rowMeans(XX)
  } else if(mode=="nodes"){
    XX<-matrix(x[tree$edge],nrow(tree$edge),2)
    x<-rowMeans(XX)
  }
  # begin optional arguments
  if(hasArg(tol)) tol<-list(...)$tol
  else tol<-1e-6
  if(hasArg(prompt)) prompt<-list(...)$prompt
  else prompt<-FALSE
  if(hasArg(type)) type<-list(...)$type
  else type<-"phylogram"
  if(hasArg(show.tip.label)) show.tip.label<-list(...)$show.tip.label
  else show.tip.label<-TRUE
  if(hasArg(show.node.label)) show.node.label<-list(...)$show.node.label
  else show.node.label<-FALSE
  if(hasArg(edge.width)) edge.width<-list(...)$edge.width
  else edge.width<-4
  if(hasArg(edge.lty)) edge.lty<-list(...)$edge.lty
  else edge.lty<-1
  if(hasArg(font)) font<-list(...)$font
  else font<-3
  if(hasArg(cex)) cex<-list(...)$cex
  else cex<-par("cex")
  if(hasArg(adj)) adj<-list(...)$adj
  else adj<-NULL
  if(hasArg(srt)) srt<-list(...)$srt
  else srt<-0
  if(hasArg(no.margin)) no.margin<-list(...)$no.margin
  else no.margin<-TRUE
  if(hasArg(label.offset)) label.offset<-list(...)$label.offset
  else label.offset<-0.01*max(nodeHeights(tree))
  if(hasArg(underscore)) underscore<-list(...)$underscore
  else underscore<-FALSE
  if(hasArg(x.lim)) x.lim<-list(...)$x.lim
  else x.lim<-NULL
  if(hasArg(y.lim)) y.lim<-list(...)$y.lim
  else y.lim<-if(legend&&!prompt) c(1-0.06*length(tree$tip.label),length(tree$tip.label)) else NULL
  if(hasArg(direction)) direction<-list(...)$direction
  else direction<-"rightwards"
  if(hasArg(lab4ut)) lab4ut<-list(...)$lab4ut
  else lab4ut<-"horizontal"
  if(hasArg(tip.color)) tip.color<-list(...)$tip.color
  else tip.color<-"black"
  if(hasArg(plot)) plot<-list(...)$plot
  else plot<-TRUE
  if(hasArg(rotate.tree)) rotate.tree<-list(...)$rotate.tree
  else rotate.tree<-0
  if(hasArg(open.angle)) open.angle<-list(...)$open.angle
  else open.angle<-0
  if(hasArg(use.edge.length)) use.edge.length<-list(...)$use.edge.length
  else use.edge.length<-TRUE
  # end optional arguments
  if(palette=="heat.colors") cols<-heat.colors(n=1000)
  if(palette=="gray") cols<-gray(1000:1/1000)
  if(palette=="rainbow")	cols<-rainbow(1000,start=0.7,end=0) # blue->red
  if(palette=="custom") {colfunct<-colorRampPalette(c("#3D52A1","#77B7E5","#FFFAD2","#ED875E","#AE1C3E"))
  cols <- colfunct(1000)}
  if(is.null(xlims)) xlims<-range(x)+c(-tol,tol)
  breaks<-0:1000/1000*(xlims[2]-xlims[1])+xlims[1]
  whichColor<-function(p,cols,breaks){
    i<-1
    while(p>=breaks[i]&&p>breaks[i+1]) i<-i+1
    cols[i]
  }
  colors<-sapply(x,whichColor,cols=cols,breaks=breaks)
  par(lend=2)
  # now plot
  tree <- compute.brlen(tree)
  xx<-plot.phylo(tree,type=type,show.tip.label=show.tip.label,show.node.label=show.node.label,edge.color=colors,
                 edge.width=edge.width,edge.lty=edge.lty,font=font,cex=cex,adj=adj,srt=srt,no.margin=no.margin,root.edge=root.edge,
                 label.offset=label.offset,underscore=underscore,x.lim=x.lim,y.lim=y.lim,direction=direction,lab4ut=lab4ut,
                 tip.color=tip.color,plot=plot,rotate.tree=rotate.tree,open.angle=open.angle,lend=2,new=FALSE,use.edge.length=FALSE)
  if(legend==TRUE&&is.logical(legend)) legend<-round(15*max(nodeHeights(tree)),2)
  if(legend){
    if(hasArg(title)) title<-list(...)$title
    else title<-"trait value"
    if(hasArg(digits)) digits<-list(...)$digits
    else digits<-1
    add.color.bar(legend,cols,title,xlims,digits,prompt=prompt)
  }
  invisible(xx)
}

sp_tree_bl_names <- sp_tree_bl
species_info_df <- species_info %>%
  as.data.frame() %>%
  column_to_rownames("code")

sp_tree_bl_names$tip.label <- species_info_df[sp_tree_bl_names$tip.label,"common_name"]

sp_names_PC1 <- pca_sp$x[,"PC1"]
names(sp_names_PC1) <- species_info_df[names(sp_names_PC1),"common_name"]
sp_names_PC2 <- pca_sp$x[,"PC2"]
names(sp_names_PC2) <- species_info_df[names(sp_names_PC2),"common_name"]
sp_names_PC3 <- pca_sp$x[,"PC3"]
names(sp_names_PC3) <- species_info_df[names(sp_names_PC3),"common_name"]

pdf("06_output_cluster_by_species/PC_gene_pval_individual_PC1.pdf")
plotBranchbyTrait(sp_tree_bl_names,sp_names_PC1,mode="tips",palette = "custom",use.edge.length=FALSE)
dev.off()
pdf("06_output_cluster_by_species/PC_gene_pval_individual_PC2.pdf")
plotBranchbyTrait(sp_tree_bl_names,sp_names_PC2,mode="tips",palette = "custom",use.edge.length=FALSE)
dev.off()
pdf("06_output_cluster_by_species/PC_gene_pval_individual_PC3.pdf")
plotBranchbyTrait(sp_tree_bl_names,sp_names_PC3,mode="tips",palette = "custom",use.edge.length=FALSE)
dev.off()




 
#########Correlations with PC axes and life history traits


sp_coord <- pca_sp$x %>%
  as.data.frame %>%
  rownames_to_column("sp_abbr") %>%
  as.tibble

sp_coord_anno <- sp_coord %>%
  left_join(species_info,by=c("sp_abbr" = "code")) %>%
  mutate(heart_index = heart_mass_mean/body_mass_mean) %>%
  mutate(log_body_mean_mean_hbabm = log(body_mass_mean_hbabm))

sp_coord_anno %>%
  ggplot(aes(log(body_mass_mean_hbabm),PC1)) +
  geom_point() +
  xlab("log(body mass)")
ggsave("06_output_cluster_by_species/PC1_gene_pval_log_BM_HABM.pdf",width=6,height=4)


sp_coord_anno_df <- sp_coord_anno %>%
  as.data.frame()
rownames(sp_coord_anno_df) <- sp_coord_anno_df$sp_abbr

#Plot phylomorphospace plot of body size by PC1
pdf("06_output_cluster_by_species/PC1_gene_pval_log_BM_HABM_phylomorphospace.pdf",width=8,height=8)
phylomorphospace(sp_tree_bl,sp_coord_anno_df[, c("log_body_mean_mean_hbabm", "PC1")], label="off", xlab="log(body mass)", ylab="PC1")
dev.off()

#Correlation between PC1 and avian body masses

bm_hbabm_pc1 <- gls(PC1 ~ log_body_mean_mean_hbabm,correlation=corBrownian(1,sp_tree_bl),data=sp_coord_anno_df)
ou_hbabm_pc1 <- gls(PC1 ~ log_body_mean_mean_hbabm,correlation=corMartins(1,sp_tree_bl),data=sp_coord_anno_df)
sink("06_output_cluster_by_species/bodymass_PC1_gene_pval_PGLSCorr.txt")
summary(bm_hbabm_pc1)
summary(ou_hbabm_pc1)
sink()


###Look at correlations with selected pvalues and body mass with Spearman's rank correlation, to account for non-normalcy (but doesn't account for phylogenetic structure)
#Is there any enrichment in KEGG pathways, given PC loading scores?
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

#Prep format for correlation
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
  as.data.frame

bm_hog_spear_res <- list()
#For each hog that has at least 1 selected branch, perform cor.test between mean body mass and selected branches. Record p-value, correlation
for (i in 1:length(names_sel_hogs)){
  hog <- names_sel_hogs[i]
  #data <- comb_data %>%
    #dplyr::select(sp_abbr,hog,body_mass_mean_hbabm) %>%
    #as.data.frame
  #rownames(data) <- data$sp_abbr
  
  x <- comb_data_df[,colnames(comb_data_df)==hog]
  y <- log(comb_data_df$body_mass_mean_hbabm)

  
  try(bm_hog_spear_res[[i]] <- cor.test(x,y,method="spearman"))
  try(hog_phylo_spear_res[i,2] <- bm_hog_spear_res[[i]]$p.value)
  try(hog_phylo_spear_res[i,3] <- bm_hog_spear_res[[i]]$estimate)
  try(hog_phylo_spear_res[i,4] <- bm_hog_spear_res[[i]]$statistic)
}

colnames(hog_phylo_spear_res) <- c("hog","pvalue","rho","statistic")

hog_phylo_spear_res <- hog_phylo_spear_res %>%
  as.tibble %>%
  mutate(pvalue = as.double(pvalue), rho = as.double(rho), statistic = as.double(statistic))

save(bm_hog_res,hog_phylo_spear_res,file ="06_output_cluster_by_species/bodymass_gene_pval_spearman_res.Rdat")
#load("06_output_cluster_by_species/bodymass_sppvals_spearman_res.Rdat")

hog_phylo_spear_res_anno <- hog_phylo_spear_res %>%
  left_join(all_res_gene_zf_hs) %>%
  dplyr::select(entrezgene,entrezgene_hs,pvalue,rho,FDRPval_busted,hog) %>%
  mutate(FDRPval_spear = p.adjust(pvalue,method="BH")) %>%
  left_join(n_sig_sp) %>%
  arrange(FDRPval_spear)



#See if any enrichment (or different PC scores) for different categories.
#Explore loadings

hog_geneids <- all_res_gene_zf_hs %>%
  dplyr::select(hog,entrezgene,entrezgene_zf,entrezgene_hs) %>%
  distinct(hog,.keep_all = TRUE) %>%
  mutate(hog = paste0("hog_",hog))


#Gene set enrichement with spearman correlation results, only test genes that are selected in at least 5 species to avoid biases associated with genes that are present in only a few species.

Corr_loading_geneList <- hog_phylo_spear_res_anno %>%
  filter(!is.na(entrezgene)) %>%
  filter(n_sel > 5) %>%
  distinct(entrezgene,.keep_all=T) %>%
  pull(rho)
names(Corr_loading_geneList) <-hog_phylo_spear_res_anno %>%
  filter(!is.na(entrezgene)) %>%
  filter(n_sel > 5) %>%
  distinct(entrezgene,.keep_all=T) %>%
  pull(entrezgene)
Corr_loading_geneList <- sort(Corr_loading_geneList,decreasing = TRUE)

corr_kk <- gseKEGG(geneList = Corr_loading_geneList,organism = "gga", keyType = "ncbi-geneid", nPerm = 1000, pvalueCutoff = 1,minGSSize = 10,by="DOSE")
corr_kk_summary <- as.data.frame(corr_kk)
write_csv(corr_kk_summary,"06_output_cluster_by_species/bodymass_PC1_gene_pval_spearman_GSE_results_nocutoff.csv")

corr_kk_summary %>%
  filter(qvalues<0.3) %>%
  write_csv("06_output_cluster_by_species/bodymass_PC1_gene_pval_spearman_GSE_results_q<0.3.csv")

#gseaplot(corr_kk,geneSetID=1,title = "cellular senescence")

#Are the results robust to looking at enrichment of genes with significant p-values? Cellular senescence results are
corr_sig <- hog_phylo_spear_res_anno %>%
  filter(FDRPval_spear < 0.3) %>%
  filter(!is.na(entrezgene)) %>%
  filter(n_sel > 5) %>%
  distinct(entrezgene) %>%
  pull(entrezgene)
all_tested <- hog_phylo_spear_res_anno %>%
  filter(n_sel > 5) %>%
  filter(!is.na(entrezgene)) %>%
  distinct(entrezgene) %>%
  pull(entrezgene)
Corr_k <- enrichKEGG(corr_sig,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=all_tested,keyType="ncbi-geneid")
spear_enrich <- as.data.frame(Corr_k)
spear_enrich[1:5,]

spear_enrich %>%
  write_csv("06_output_cluster_by_species/bodymass_PC1_gene_pval_spearman_enrich_results_nocutoff.csv")

spear_enrich %>%
  filter(qvalue<0.2) %>%
  write_csv("06_output_cluster_by_species/bodymass_PC1_gene_pval_spearman_enrich_results_q<0.2.csv")






##############Do all analyses using omega

pca_sp <- pca_omega

#Read in tree, chose random hog to give branch lengths. Species tree has same topology for all hogs. Note that results are robust to hog selected.
sp_tree_bl <- read.tree("06_input_cluster_by_species/11700.tree1.nwk")
#sp_tree_bl <- read.tree("Alignments/10090.tree1.nwk")
sp_tree <- read.tree("06_input_cluster_by_species/11700.final_spt.nwk")

#Making custom version of plotBranchbyTrait from phytools to include plotting without branch lengths and custom color palette.
plotBranchbyTrait<-function(tree,x,mode=c("edges","tips","nodes"),palette="rainbow",legend=TRUE,xlims=NULL,...){
  mode<-mode[1]
  if(mode=="tips"){
    x<-c(x[tree$tip.label],fastAnc(tree,x))
    names(x)[1:length(tree$tip.label)]<-1:length(tree$tip.label)
    XX<-matrix(x[tree$edge],nrow(tree$edge),2)
    x<-rowMeans(XX)
  } else if(mode=="nodes"){
    XX<-matrix(x[tree$edge],nrow(tree$edge),2)
    x<-rowMeans(XX)
  }
  # begin optional arguments
  if(hasArg(tol)) tol<-list(...)$tol
  else tol<-1e-6
  if(hasArg(prompt)) prompt<-list(...)$prompt
  else prompt<-FALSE
  if(hasArg(type)) type<-list(...)$type
  else type<-"phylogram"
  if(hasArg(show.tip.label)) show.tip.label<-list(...)$show.tip.label
  else show.tip.label<-TRUE
  if(hasArg(show.node.label)) show.node.label<-list(...)$show.node.label
  else show.node.label<-FALSE
  if(hasArg(edge.width)) edge.width<-list(...)$edge.width
  else edge.width<-4
  if(hasArg(edge.lty)) edge.lty<-list(...)$edge.lty
  else edge.lty<-1
  if(hasArg(font)) font<-list(...)$font
  else font<-3
  if(hasArg(cex)) cex<-list(...)$cex
  else cex<-par("cex")
  if(hasArg(adj)) adj<-list(...)$adj
  else adj<-NULL
  if(hasArg(srt)) srt<-list(...)$srt
  else srt<-0
  if(hasArg(no.margin)) no.margin<-list(...)$no.margin
  else no.margin<-TRUE
  if(hasArg(label.offset)) label.offset<-list(...)$label.offset
  else label.offset<-0.01*max(nodeHeights(tree))
  if(hasArg(underscore)) underscore<-list(...)$underscore
  else underscore<-FALSE
  if(hasArg(x.lim)) x.lim<-list(...)$x.lim
  else x.lim<-NULL
  if(hasArg(y.lim)) y.lim<-list(...)$y.lim
  else y.lim<-if(legend&&!prompt) c(1-0.06*length(tree$tip.label),length(tree$tip.label)) else NULL
  if(hasArg(direction)) direction<-list(...)$direction
  else direction<-"rightwards"
  if(hasArg(lab4ut)) lab4ut<-list(...)$lab4ut
  else lab4ut<-"horizontal"
  if(hasArg(tip.color)) tip.color<-list(...)$tip.color
  else tip.color<-"black"
  if(hasArg(plot)) plot<-list(...)$plot
  else plot<-TRUE
  if(hasArg(rotate.tree)) rotate.tree<-list(...)$rotate.tree
  else rotate.tree<-0
  if(hasArg(open.angle)) open.angle<-list(...)$open.angle
  else open.angle<-0
  if(hasArg(use.edge.length)) use.edge.length<-list(...)$use.edge.length
  else use.edge.length<-TRUE
  # end optional arguments
  if(palette=="heat.colors") cols<-heat.colors(n=1000)
  if(palette=="gray") cols<-gray(1000:1/1000)
  if(palette=="rainbow")	cols<-rainbow(1000,start=0.7,end=0) # blue->red
  if(palette=="custom") {colfunct<-colorRampPalette(c("#3D52A1","#77B7E5","#FFFAD2","#ED875E","#AE1C3E"))
  cols <- colfunct(1000)}
  if(is.null(xlims)) xlims<-range(x)+c(-tol,tol)
  breaks<-0:1000/1000*(xlims[2]-xlims[1])+xlims[1]
  whichColor<-function(p,cols,breaks){
    i<-1
    while(p>=breaks[i]&&p>breaks[i+1]) i<-i+1
    cols[i]
  }
  colors<-sapply(x,whichColor,cols=cols,breaks=breaks)
  par(lend=2)
  # now plot
  tree <- compute.brlen(tree)
  xx<-plot.phylo(tree,type=type,show.tip.label=show.tip.label,show.node.label=show.node.label,edge.color=colors,
                 edge.width=edge.width,edge.lty=edge.lty,font=font,cex=cex,adj=adj,srt=srt,no.margin=no.margin,root.edge=root.edge,
                 label.offset=label.offset,underscore=underscore,x.lim=x.lim,y.lim=y.lim,direction=direction,lab4ut=lab4ut,
                 tip.color=tip.color,plot=plot,rotate.tree=rotate.tree,open.angle=open.angle,lend=2,new=FALSE,use.edge.length=FALSE)
  if(legend==TRUE&&is.logical(legend)) legend<-round(15*max(nodeHeights(tree)),2)
  if(legend){
    if(hasArg(title)) title<-list(...)$title
    else title<-"trait value"
    if(hasArg(digits)) digits<-list(...)$digits
    else digits<-1
    add.color.bar(legend,cols,title,xlims,digits,prompt=prompt)
  }
  invisible(xx)
}

sp_tree_bl_names <- sp_tree_bl
species_info_df <- species_info %>%
  as.data.frame() %>%
  column_to_rownames("code")

sp_tree_bl_names$tip.label <- species_info_df[sp_tree_bl_names$tip.label,"common_name"]

sp_names_PC1 <- pca_sp$x[,"PC1"]
names(sp_names_PC1) <- species_info_df[names(sp_names_PC1),"common_name"]
sp_names_PC2 <- pca_sp$x[,"PC2"]
names(sp_names_PC2) <- species_info_df[names(sp_names_PC2),"common_name"]
sp_names_PC3 <- pca_sp$x[,"PC3"]
names(sp_names_PC3) <- species_info_df[names(sp_names_PC3),"common_name"]

pdf("06_output_cluster_by_species/PC_gene_omega_individual_PC1.pdf")
plotBranchbyTrait(sp_tree_bl_names,sp_names_PC1,mode="tips",palette = "custom",use.edge.length=FALSE)
dev.off()
pdf("06_output_cluster_by_species/PC_gene_omega_individual_PC2.pdf")
plotBranchbyTrait(sp_tree_bl_names,sp_names_PC2,mode="tips",palette = "custom",use.edge.length=FALSE)
dev.off()
pdf("06_output_cluster_by_species/PC_gene_omega_individual_PC3.pdf")
plotBranchbyTrait(sp_tree_bl_names,sp_names_PC3,mode="tips",palette = "custom",use.edge.length=FALSE)
dev.off()





#########Correlations with PC axes and life history traits


sp_coord <- pca_sp$x %>%
  as.data.frame %>%
  rownames_to_column("sp_abbr") %>%
  as.tibble

sp_coord_anno <- sp_coord %>%
  left_join(species_info,by=c("sp_abbr" = "code")) %>%
  mutate(heart_index = heart_mass_mean/body_mass_mean) %>%
  mutate(log_body_mean_mean_hbabm = log(body_mass_mean_hbabm))

sp_coord_anno %>%
  ggplot(aes(log(body_mass_mean_hbabm),PC1)) +
  geom_point() +
  xlab("log(body mass)")
ggsave("06_output_cluster_by_species/PC1_gene_omega_log_BM_HABM.pdf",width=6,height=4)


sp_coord_anno_df <- sp_coord_anno %>%
  as.data.frame()
rownames(sp_coord_anno_df) <- sp_coord_anno_df$sp_abbr

#Plot phylomorphospace plot of body size by PC1
pdf("06_output_cluster_by_species/PC1_gene_omega_log_BM_HABM_phylomorphospace.pdf",width=8,height=8)
phylomorphospace(sp_tree_bl,sp_coord_anno_df[, c("log_body_mean_mean_hbabm", "PC1")], label="off", xlab="log(body mass)", ylab="PC1")
dev.off()


#Correlation between PC1 and avian body masses

bm_hbabm_pc1 <- gls(PC1 ~ log_body_mean_mean_hbabm,correlation=corBrownian(1,sp_tree_bl),data=sp_coord_anno_df)
ou_hbabm_pc1 <- gls(PC1 ~ log_body_mean_mean_hbabm,correlation=corMartins(1,sp_tree_bl),data=sp_coord_anno_df)
sink("06_output_cluster_by_species/bodymass_PC1_gene_omega_PGLSCorr.txt")
summary(bm_hbabm_pc1)
summary(ou_hbabm_pc1)
sink()


###Look at correlations with selected pvalues and body mass with Spearman's rank correlation, to account for non-normalcy (but doesn't account for phylogenetic structure)
#Is there any enrichment in KEGG pathways, given PC loading scores?
bsrel_sp_pval_res_gene_raw_na <- bsrel_sp_pval_res_gene %>%
  filter(pval_type == "pval") %>%
  separate(hog_treenum,into=c("hog","treenum")) %>%
  dplyr::select(-pval_type,-treenum) %>%
  gather(sp,pval,-hog) %>%
  mutate(pval=ifelse(pval==0,1e-18,pval)) %>%
  mutate(pval=log10(pval)) %>%
  mutate(hog=paste0("hog_",hog)) %>%
  spread(hog,pval)

bsrel_sp_omega_res_gene_raw_na <- bsrel_sp_all_params_res_gene_raw %>%
  dplyr::select(hog,sp,omega_sig) %>%
  mutate(omega_sig=ifelse(omega_sig==0,1,omega_sig)) %>%
  mutate(omega_sig=log10(omega_sig)) %>%
  spread(hog,omega_sig)

bsrel_sp_weights_res_gene_raw_na <- bsrel_sp_all_params_res_gene_raw %>%
  dplyr::select(hog,sp,weight_sig) %>%
  spread(hog,weight_sig)

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

#Prep format for correlation omega
bsrel_sp_omega_res_gene_forcorr <- bsrel_sp_omega_res_gene_raw_na %>%
  dplyr::select(sp,hogs_to_keep) %>%
  as.data.frame
rownames(bsrel_sp_omega_res_gene_forcorr) <- bsrel_sp_omega_res_gene_forcorr$sp
bsrel_sp_omega_res_gene_forcorr$sp <- NULL



#Pull vector of hogs that have at least 1 selected branch:
sel_hogs <- n_sig_sp %>% filter(n_sel > 1) %>% pull(hog)
names_sel_hogs <- paste0("hog_",sel_hogs)

hog_phylo_spear_res <- matrix(nrow=length(names_sel_hogs),ncol=4)
hog_phylo_spear_res[,1] <- sel_hogs


#Create a combo dataset with selected branches and species data
comb_data <- bsrel_sp_omega_res_gene_forcorr %>%
  rownames_to_column(var = "sp_abbr") %>%
  as.tibble %>%
  left_join(sp_coord_anno)

comb_data_df <- comb_data %>%
  as.data.frame

bm_hog_spear_res <- list()
#For each hog that has at least 1 selected branch, perform cor.test between mean body mass and selected branches. Record p-value, correlation
for (i in 1:length(names_sel_hogs)){
  hog <- names_sel_hogs[i]
  #data <- comb_data %>%
  #dplyr::select(sp_abbr,hog,body_mass_mean_hbabm) %>%
  #as.data.frame
  #rownames(data) <- data$sp_abbr
  
  x <- comb_data_df[,colnames(comb_data_df)==hog]
  y <- log(comb_data_df$body_mass_mean_hbabm)
  
  
  try(bm_hog_spear_res[[i]] <- cor.test(x,y,method="spearman"))
  try(hog_phylo_spear_res[i,2] <- bm_hog_spear_res[[i]]$p.value)
  try(hog_phylo_spear_res[i,3] <- bm_hog_spear_res[[i]]$estimate)
  try(hog_phylo_spear_res[i,4] <- bm_hog_spear_res[[i]]$statistic)
}

colnames(hog_phylo_spear_res) <- c("hog","pvalue","rho","statistic")

hog_phylo_spear_res <- hog_phylo_spear_res %>%
  as.tibble %>%
  mutate(pvalue = as.double(pvalue), rho = as.double(rho), statistic = as.double(statistic))

save(bm_hog_res,hog_phylo_spear_res,file ="06_output_cluster_by_species/bodymass_gene_omega_spearman_res.Rdat")
#load("06_output_cluster_by_species/bodymass_sppvals_spearman_res.Rdat")

hog_phylo_spear_res_anno <- hog_phylo_spear_res %>%
  left_join(all_res_gene_zf_hs) %>%
  dplyr::select(entrezgene,entrezgene_hs,pvalue,rho,FDRPval_busted,hog) %>%
  mutate(FDRPval_spear = p.adjust(pvalue,method="BH")) %>%
  left_join(n_sig_sp) %>%
  arrange(FDRPval_spear)



#See if any enrichment (or different PC scores) for different categories.
#Explore loadings

hog_geneids <- all_res_gene_zf_hs %>%
  dplyr::select(hog,entrezgene,entrezgene_zf,entrezgene_hs) %>%
  distinct(hog,.keep_all = TRUE) %>%
  mutate(hog = paste0("hog_",hog))


#Gene set enrichement with spearman correlation results, only test genes that are selected in at least 5 species to avoid biases associated with genes that are present in only a few species.

Corr_loading_geneList <- hog_phylo_spear_res_anno %>%
  filter(!is.na(entrezgene)) %>%
  filter(n_sel > 5) %>%
  distinct(entrezgene,.keep_all=T) %>%
  pull(rho)
names(Corr_loading_geneList) <-hog_phylo_spear_res_anno %>%
  filter(!is.na(entrezgene)) %>%
  filter(n_sel > 5) %>%
  distinct(entrezgene,.keep_all=T) %>%
  pull(entrezgene)
Corr_loading_geneList <- sort(Corr_loading_geneList,decreasing = TRUE)

corr_kk <- gseKEGG(geneList = Corr_loading_geneList,organism = "gga", keyType = "ncbi-geneid", nPerm = 1000, pvalueCutoff = 1,minGSSize = 10,by="DOSE")
corr_kk_summary <- as.data.frame(corr_kk)
write_csv(corr_kk_summary,"06_output_cluster_by_species/bodymass_PC1_gene_omega_spearman_GSE_results_nocutoff.csv")

corr_kk_summary %>%
  filter(qvalues<0.3) %>%
  write_csv("06_output_cluster_by_species/bodymass_PC1_gene_omega_spearman_GSE_results_q<0.3.csv")

#gseaplot(corr_kk,geneSetID=1,title = "cellular senescence")

#Are the results robust to looking at enrichment of genes with significant p-values? Cellular senescence results are
corr_sig <- hog_phylo_spear_res_anno %>%
  filter(FDRPval_spear < 0.3) %>%
  filter(!is.na(entrezgene)) %>%
  filter(n_sel > 5) %>%
  distinct(entrezgene) %>%
  pull(entrezgene)
all_tested <- hog_phylo_spear_res_anno %>%
  filter(n_sel > 5) %>%
  filter(!is.na(entrezgene)) %>%
  distinct(entrezgene) %>%
  pull(entrezgene)
Corr_k <- enrichKEGG(corr_sig,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=all_tested,keyType="ncbi-geneid")
spear_enrich <- as.data.frame(Corr_k)
spear_enrich[1:5,]

spear_enrich %>%
  write_csv("06_output_cluster_by_species/bodymass_PC1_gene_omega_spearman_enrich_results_nocutoff.csv")

spear_enrich %>%
  filter(qvalue<0.2) %>%
  write_csv("06_output_cluster_by_species/bodymass_PC1_gene_omega_spearman_enrich_results_q<0.2.csv")








###########Do it all using weights
pca_sp <- pca_weights

#Read in tree, chose random hog to give branch lengths. Species tree has same topology for all hogs. Note that results are robust to hog selected.
sp_tree_bl <- read.tree("06_input_cluster_by_species/11700.tree1.nwk")
#sp_tree_bl <- read.tree("Alignments/10090.tree1.nwk")
sp_tree <- read.tree("06_input_cluster_by_species/11700.final_spt.nwk")

#Making custom version of plotBranchbyTrait from phytools to include plotting without branch lengths and custom color palette.
plotBranchbyTrait<-function(tree,x,mode=c("edges","tips","nodes"),palette="rainbow",legend=TRUE,xlims=NULL,...){
  mode<-mode[1]
  if(mode=="tips"){
    x<-c(x[tree$tip.label],fastAnc(tree,x))
    names(x)[1:length(tree$tip.label)]<-1:length(tree$tip.label)
    XX<-matrix(x[tree$edge],nrow(tree$edge),2)
    x<-rowMeans(XX)
  } else if(mode=="nodes"){
    XX<-matrix(x[tree$edge],nrow(tree$edge),2)
    x<-rowMeans(XX)
  }
  # begin optional arguments
  if(hasArg(tol)) tol<-list(...)$tol
  else tol<-1e-6
  if(hasArg(prompt)) prompt<-list(...)$prompt
  else prompt<-FALSE
  if(hasArg(type)) type<-list(...)$type
  else type<-"phylogram"
  if(hasArg(show.tip.label)) show.tip.label<-list(...)$show.tip.label
  else show.tip.label<-TRUE
  if(hasArg(show.node.label)) show.node.label<-list(...)$show.node.label
  else show.node.label<-FALSE
  if(hasArg(edge.width)) edge.width<-list(...)$edge.width
  else edge.width<-4
  if(hasArg(edge.lty)) edge.lty<-list(...)$edge.lty
  else edge.lty<-1
  if(hasArg(font)) font<-list(...)$font
  else font<-3
  if(hasArg(cex)) cex<-list(...)$cex
  else cex<-par("cex")
  if(hasArg(adj)) adj<-list(...)$adj
  else adj<-NULL
  if(hasArg(srt)) srt<-list(...)$srt
  else srt<-0
  if(hasArg(no.margin)) no.margin<-list(...)$no.margin
  else no.margin<-TRUE
  if(hasArg(label.offset)) label.offset<-list(...)$label.offset
  else label.offset<-0.01*max(nodeHeights(tree))
  if(hasArg(underscore)) underscore<-list(...)$underscore
  else underscore<-FALSE
  if(hasArg(x.lim)) x.lim<-list(...)$x.lim
  else x.lim<-NULL
  if(hasArg(y.lim)) y.lim<-list(...)$y.lim
  else y.lim<-if(legend&&!prompt) c(1-0.06*length(tree$tip.label),length(tree$tip.label)) else NULL
  if(hasArg(direction)) direction<-list(...)$direction
  else direction<-"rightwards"
  if(hasArg(lab4ut)) lab4ut<-list(...)$lab4ut
  else lab4ut<-"horizontal"
  if(hasArg(tip.color)) tip.color<-list(...)$tip.color
  else tip.color<-"black"
  if(hasArg(plot)) plot<-list(...)$plot
  else plot<-TRUE
  if(hasArg(rotate.tree)) rotate.tree<-list(...)$rotate.tree
  else rotate.tree<-0
  if(hasArg(open.angle)) open.angle<-list(...)$open.angle
  else open.angle<-0
  if(hasArg(use.edge.length)) use.edge.length<-list(...)$use.edge.length
  else use.edge.length<-TRUE
  # end optional arguments
  if(palette=="heat.colors") cols<-heat.colors(n=1000)
  if(palette=="gray") cols<-gray(1000:1/1000)
  if(palette=="rainbow")	cols<-rainbow(1000,start=0.7,end=0) # blue->red
  if(palette=="custom") {colfunct<-colorRampPalette(c("#3D52A1","#77B7E5","#FFFAD2","#ED875E","#AE1C3E"))
  cols <- colfunct(1000)}
  if(is.null(xlims)) xlims<-range(x)+c(-tol,tol)
  breaks<-0:1000/1000*(xlims[2]-xlims[1])+xlims[1]
  whichColor<-function(p,cols,breaks){
    i<-1
    while(p>=breaks[i]&&p>breaks[i+1]) i<-i+1
    cols[i]
  }
  colors<-sapply(x,whichColor,cols=cols,breaks=breaks)
  par(lend=2)
  # now plot
  tree <- compute.brlen(tree)
  xx<-plot.phylo(tree,type=type,show.tip.label=show.tip.label,show.node.label=show.node.label,edge.color=colors,
                 edge.width=edge.width,edge.lty=edge.lty,font=font,cex=cex,adj=adj,srt=srt,no.margin=no.margin,root.edge=root.edge,
                 label.offset=label.offset,underscore=underscore,x.lim=x.lim,y.lim=y.lim,direction=direction,lab4ut=lab4ut,
                 tip.color=tip.color,plot=plot,rotate.tree=rotate.tree,open.angle=open.angle,lend=2,new=FALSE,use.edge.length=FALSE)
  if(legend==TRUE&&is.logical(legend)) legend<-round(15*max(nodeHeights(tree)),2)
  if(legend){
    if(hasArg(title)) title<-list(...)$title
    else title<-"trait value"
    if(hasArg(digits)) digits<-list(...)$digits
    else digits<-1
    add.color.bar(legend,cols,title,xlims,digits,prompt=prompt)
  }
  invisible(xx)
}

sp_tree_bl_names <- sp_tree_bl
species_info_df <- species_info %>%
  as.data.frame() %>%
  column_to_rownames("code")

sp_tree_bl_names$tip.label <- species_info_df[sp_tree_bl_names$tip.label,"common_name"]

sp_names_PC1 <- pca_sp$x[,"PC1"]
names(sp_names_PC1) <- species_info_df[names(sp_names_PC1),"common_name"]
sp_names_PC2 <- pca_sp$x[,"PC2"]
names(sp_names_PC2) <- species_info_df[names(sp_names_PC2),"common_name"]
sp_names_PC3 <- pca_sp$x[,"PC3"]
names(sp_names_PC3) <- species_info_df[names(sp_names_PC3),"common_name"]

pdf("06_output_cluster_by_species/PC_gene_weights_individual_PC1.pdf")
plotBranchbyTrait(sp_tree_bl_names,sp_names_PC1,mode="tips",palette = "custom",use.edge.length=FALSE)
dev.off()
pdf("06_output_cluster_by_species/PC_gene_weights_individual_PC2.pdf")
plotBranchbyTrait(sp_tree_bl_names,sp_names_PC2,mode="tips",palette = "custom",use.edge.length=FALSE)
dev.off()
pdf("06_output_cluster_by_species/PC_gene_weights_individual_PC3.pdf")
plotBranchbyTrait(sp_tree_bl_names,sp_names_PC3,mode="tips",palette = "custom",use.edge.length=FALSE)
dev.off()





#########Correlations with PC axes and life history traits


sp_coord <- pca_sp$x %>%
  as.data.frame %>%
  rownames_to_column("sp_abbr") %>%
  as.tibble

sp_coord_anno <- sp_coord %>%
  left_join(species_info,by=c("sp_abbr" = "code")) %>%
  mutate(heart_index = heart_mass_mean/body_mass_mean) %>%
  mutate(log_body_mean_mean_hbabm = log(body_mass_mean_hbabm))

sp_coord_anno %>%
  ggplot(aes(log(body_mass_mean_hbabm),PC1)) +
  geom_point() +
  xlab("log(body mass)")
ggsave("06_output_cluster_by_species/PC1_gene_weights_log_BM_HABM.pdf",width=6,height=4)


sp_coord_anno_df <- sp_coord_anno %>%
  as.data.frame()
rownames(sp_coord_anno_df) <- sp_coord_anno_df$sp_abbr

#Plot phylomorphospace plot of body size by PC1
pdf("06_output_cluster_by_species/PC1_gene_weightslog_BM_HABM_phylomorphospace.pdf",width=8,height=8)
phylomorphospace(sp_tree_bl,sp_coord_anno_df[, c("log_body_mean_mean_hbabm", "PC1")], label="off", xlab="log(body mass)", ylab="PC1")
dev.off()



#Correlation between PC1 and avian body masses

bm_hbabm_pc1 <- gls(PC1 ~ log_body_mean_mean_hbabm,correlation=corBrownian(1,sp_tree_bl),data=sp_coord_anno_df)
ou_hbabm_pc1 <- gls(PC1 ~ log_body_mean_mean_hbabm,correlation=corMartins(1,sp_tree_bl),data=sp_coord_anno_df)
sink("06_output_cluster_by_species/bodymass_PC1_gene_weights_PGLSCorr.txt")
summary(bm_hbabm_pc1)
summary(ou_hbabm_pc1)
sink()


###Look at correlations with selected pvalues and body mass with Spearman's rank correlation, to account for non-normalcy (but doesn't account for phylogenetic structure)
#Is there any enrichment in KEGG pathways, given PC loading scores?
bsrel_sp_pval_res_gene_raw_na <- bsrel_sp_pval_res_gene %>%
  filter(pval_type == "pval") %>%
  separate(hog_treenum,into=c("hog","treenum")) %>%
  dplyr::select(-pval_type,-treenum) %>%
  gather(sp,pval,-hog) %>%
  mutate(pval=ifelse(pval==0,1e-18,pval)) %>%
  mutate(pval=log10(pval)) %>%
  mutate(hog=paste0("hog_",hog)) %>%
  spread(hog,pval)

bsrel_sp_omega_res_gene_raw_na <- bsrel_sp_all_params_res_gene_raw %>%
  dplyr::select(hog,sp,omega_sig) %>%
  mutate(omega_sig=ifelse(omega_sig==0,1,omega_sig)) %>%
  mutate(omega_sig=log10(omega_sig)) %>%
  spread(hog,omega_sig)

bsrel_sp_weights_res_gene_raw_na <- bsrel_sp_all_params_res_gene_raw %>%
  dplyr::select(hog,sp,weight_sig) %>%
  spread(hog,weight_sig)

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


#Prep format for correlation weights
bsrel_sp_weights_res_gene_forcorr <- bsrel_sp_weights_res_gene_raw_na %>%
  dplyr::select(sp,hogs_to_keep) %>%
  as.data.frame
rownames(bsrel_sp_weights_res_gene_forcorr) <- bsrel_sp_weights_res_gene_forcorr$sp
bsrel_sp_weights_res_gene_forcorr$sp <- NULL


#Pull vector of hogs that have at least 1 selected branch:
sel_hogs <- n_sig_sp %>% filter(n_sel > 1) %>% pull(hog)
names_sel_hogs <- paste0("hog_",sel_hogs)

hog_phylo_spear_res <- matrix(nrow=length(names_sel_hogs),ncol=4)
hog_phylo_spear_res[,1] <- sel_hogs

#Create a combo dataset with selected branches and species data
comb_data <- bsrel_sp_weights_res_gene_forcorr %>%
  rownames_to_column(var = "sp_abbr") %>%
  as.tibble %>%
  left_join(sp_coord_anno)

comb_data_df <- comb_data %>%
  as.data.frame

bm_hog_spear_res <- list()
#For each hog that has at least 1 selected branch, perform cor.test between mean body mass and selected branches. Record p-value, correlation
for (i in 1:length(names_sel_hogs)){
  hog <- names_sel_hogs[i]
  #data <- comb_data %>%
  #dplyr::select(sp_abbr,hog,body_mass_mean_hbabm) %>%
  #as.data.frame
  #rownames(data) <- data$sp_abbr
  
  x <- comb_data_df[,colnames(comb_data_df)==hog]
  y <- log(comb_data_df$body_mass_mean_hbabm)
  
  
  try(bm_hog_spear_res[[i]] <- cor.test(x,y,method="spearman"))
  try(hog_phylo_spear_res[i,2] <- bm_hog_spear_res[[i]]$p.value)
  try(hog_phylo_spear_res[i,3] <- bm_hog_spear_res[[i]]$estimate)
  try(hog_phylo_spear_res[i,4] <- bm_hog_spear_res[[i]]$statistic)
}

colnames(hog_phylo_spear_res) <- c("hog","pvalue","rho","statistic")

hog_phylo_spear_res <- hog_phylo_spear_res %>%
  as.tibble %>%
  mutate(pvalue = as.double(pvalue), rho = as.double(rho), statistic = as.double(statistic))

save(bm_hog_res,hog_phylo_spear_res,file ="06_output_cluster_by_species/bodymass_gene_weights_spearman_res.Rdat")
#load("06_output_cluster_by_species/bodymass_sppvals_spearman_res.Rdat")

hog_phylo_spear_res_anno <- hog_phylo_spear_res %>%
  left_join(all_res_gene_zf_hs) %>%
  dplyr::select(entrezgene,entrezgene_hs,pvalue,rho,FDRPval_busted,hog) %>%
  mutate(FDRPval_spear = p.adjust(pvalue,method="BH")) %>%
  left_join(n_sig_sp) %>%
  arrange(FDRPval_spear)



#See if any enrichment (or different PC scores) for different categories.
#Explore loadings

hog_geneids <- all_res_gene_zf_hs %>%
  dplyr::select(hog,entrezgene,entrezgene_zf,entrezgene_hs) %>%
  distinct(hog,.keep_all = TRUE) %>%
  mutate(hog = paste0("hog_",hog))


#Gene set enrichement with spearman correlation results, only test genes that are selected in at least 5 species to avoid biases associated with genes that are present in only a few species.

Corr_loading_geneList <- hog_phylo_spear_res_anno %>%
  filter(!is.na(entrezgene)) %>%
  filter(n_sel > 5) %>%
  distinct(entrezgene,.keep_all=T) %>%
  pull(rho)
names(Corr_loading_geneList) <-hog_phylo_spear_res_anno %>%
  filter(!is.na(entrezgene)) %>%
  filter(n_sel > 5) %>%
  distinct(entrezgene,.keep_all=T) %>%
  pull(entrezgene)
Corr_loading_geneList <- sort(Corr_loading_geneList,decreasing = TRUE)

corr_kk <- gseKEGG(geneList = Corr_loading_geneList,organism = "gga", keyType = "ncbi-geneid", nPerm = 1000, pvalueCutoff = 1,minGSSize = 10,by="DOSE")
corr_kk_summary <- as.data.frame(corr_kk)
write_csv(corr_kk_summary,"06_output_cluster_by_species/bodymass_PC1_gene_weights_spearman_GSE_results_nocutoff.csv")

corr_kk_summary %>%
  filter(qvalues<0.3) %>%
  write_csv("06_output_cluster_by_species/bodymass_PC1_gene_weights_spearman_GSE_results_q<0.3.csv")

#gseaplot(corr_kk,geneSetID=1,title = "cellular senescence")

#Are the results robust to looking at enrichment of genes with significant p-values? Cellular senescence results are
corr_sig <- hog_phylo_spear_res_anno %>%
  filter(FDRPval_spear < 0.3) %>%
  filter(!is.na(entrezgene)) %>%
  filter(n_sel > 5) %>%
  distinct(entrezgene) %>%
  pull(entrezgene)
all_tested <- hog_phylo_spear_res_anno %>%
  filter(n_sel > 5) %>%
  filter(!is.na(entrezgene)) %>%
  distinct(entrezgene) %>%
  pull(entrezgene)
Corr_k <- enrichKEGG(corr_sig,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=all_tested,keyType="ncbi-geneid")
spear_enrich <- as.data.frame(Corr_k)
spear_enrich[1:5,]

spear_enrich %>%
  write_csv("06_output_cluster_by_species/bodymass_PC1_gene_weights_spearman_enrich_results_nocutoff.csv")

spear_enrich %>%
  filter(qvalue<0.2) %>%
  write_csv("06_output_cluster_by_species/bodymass_PC1_gene_weights_spearman_enrich_results_q<0.2.csv")
