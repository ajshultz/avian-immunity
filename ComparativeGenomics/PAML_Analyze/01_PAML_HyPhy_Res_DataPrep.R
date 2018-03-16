setwd("~/Dropbox/BirdImmuneGeneEvolution")
library(tidyverse)

#Read in results from original analyses, and all hogs that had to be rerun, change model_num to character and correct 2a or 8a model when appropriate.
#Tree1
paml_res <- read_tsv("01_input_raw_paml_hyphy_results/paml_res_allhogs_rewrite.txt") %>% mutate(model_num = as.character(model_num))
paml_res_rerun <- read_tsv("01_input_raw_paml_hyphy_results/hogs_to_rerun_results.txt") %>% mutate(model_num = as.character(model_num))
paml_res_rerun_m7 <- read_tsv("01_input_raw_paml_hyphy_results/hogs_to_rerun_results_m7Only.txt") %>% mutate(model_num = as.character(model_num))
paml_res_rerun_m8 <- read_tsv("01_input_raw_paml_hyphy_results/hogs_to_rerun_results_m8Only.txt") %>% mutate(model_num = as.character(model_num))
paml_res_rerun_m1 <- read_tsv("01_input_raw_paml_hyphy_results/m1_rerun_hogs_results.txt") %>% mutate(model_num = as.character(model_num))
paml_res_rerun_m2 <- read_tsv("01_input_raw_paml_hyphy_results/m2_rerun_hogs_results.txt") %>% mutate(model_num = as.character(model_num))
paml_res_rerun_m7_2 <- read_tsv("01_input_raw_paml_hyphy_results/m7_rerun_hogs_results.txt") %>% mutate(model_num = as.character(model_num))
paml_res_rerun_m7_3 <- read_tsv("01_input_raw_paml_hyphy_results/m7_rerun_hogs_results_2.txt") %>% mutate(model_num = as.character(model_num))
paml_res_rerun_m8_2 <- read_tsv("01_input_raw_paml_hyphy_results/m8_rerun_hogs_results.txt") %>% mutate(model_num = as.character(model_num))
paml_res_rerun_m8_3 <- read_tsv("01_input_raw_paml_hyphy_results/m8_rerun_hogs_results_2.txt") %>% mutate(model_num = as.character(model_num))
paml_res_m2a <- read_tsv("01_input_raw_paml_hyphy_results/allhogs_results_m2aOnly.txt") %>% mutate(model_num = "2a")
paml_res_m8a <- read_tsv("01_input_raw_paml_hyphy_results/allhogs_results_m8aOnly.txt") %>% mutate(model_num = "8a")
paml_res_rerun_m2a <- read_tsv("01_input_raw_paml_hyphy_results/m2a_rerun_hogs_results.txt") %>% mutate(model_num = "2a")
paml_res_rerun_m2a_2 <- read_tsv("01_input_raw_paml_hyphy_results/m2a_rerun_hogs_results_2.txt") %>% mutate(model_num = "2a")
paml_res_rerun_m8a <- read_tsv("01_input_raw_paml_hyphy_results/m8a_rerun_hogs_results.txt") %>% mutate(model_num = "8a")
paml_res_rerun_m8a_2 <- read_tsv("01_input_raw_paml_hyphy_results/m8a_rerun_hogs_results_2.txt") %>% mutate(model_num = "8a")
paml_res_rerun_m8a_3 <- read_tsv("01_input_raw_paml_hyphy_results/m8a_rerun_hogs_results_3.txt") %>% mutate(model_num = "8a")

#Tree2
paml_res_t2 <- read_tsv("01_input_raw_paml_hyphy_results/hogs_t2_m1278_results.txt") %>% mutate(model_num = as.character(model_num))
paml_res_t2_m2a <- read_tsv("01_input_raw_paml_hyphy_results/hogs_t2_m2a_results.txt") %>% mutate(model_num = "2a")
paml_res_t2_m8a <- read_tsv("01_input_raw_paml_hyphy_results/hogs_t2_m8a_results.txt") %>% mutate(model_num = "8a")
paml_res_t2_m7_2 <- read_tsv("01_input_raw_paml_hyphy_results/hogs_t2_m7_2.txt") %>% mutate(model_num = as.character(model_num))
paml_res_t2_m8_2 <- read_tsv("01_input_raw_paml_hyphy_results/hogs_t2_m8_2.txt") %>% mutate(model_num = as.character(model_num))
paml_res_t2_m8a_2 <- read_tsv("01_input_raw_paml_hyphy_results/hogs_t2_m8a_2.txt") %>% mutate(model_num = "8a")

#Read all m0 data
paml_res_m0 <- read.table("01_input_raw_paml_hyphy_results/paml_M0_parsed.txt",header=T)
colnames(paml_res_m0) <- c("hog","model","treenum","species_tree","newick_string","lnl","treelen","kappa_m0","omega_m0")
paml_res_m0 <- paml_res_m0 %>%
  tbl_df %>%
  unite(hog,treenum,col=hog_treenum) %>%
  dplyr::select(hog_treenum,species_tree,newick_string,lnl_m0=lnl,treelen_m0=treelen,kappa_m0,omega_m0)

#Combine all PAML results:
paml_res_t1 <- paml_res %>% bind_rows(paml_res_rerun,paml_res_rerun_m7,paml_res_rerun_m8,paml_res_rerun_m1,paml_res_rerun_m2,paml_res_rerun_m7_2,paml_res_rerun_m8_2,paml_res_m2a,paml_res_m8a,paml_res_rerun_m2a,paml_res_rerun_m8a,paml_res_rerun_m8a_2,paml_res_rerun_m2a_2,paml_res_rerun_m7_3,paml_res_rerun_m8_3,paml_res_rerun_m8a_3)

paml_res_t2 <- paml_res_t2 %>% bind_rows(paml_res_t2_m2a,paml_res_t2_m8a,paml_res_t2_m7_2,paml_res_t2_m8_2,paml_res_t2_m8a_2)

#Add both tree1 and tree2 data together, add m0 results
all_paml <- paml_res_t1 %>%
  bind_rows(paml_res_t2) %>%
  unite(hog,treenum,col=hog_treenum,remove=FALSE) %>%
  left_join(paml_res_m0,by="hog_treenum") %>%
  dplyr::select(-species_tree.x,-newick_string.x) %>%
  dplyr::rename(species_tree=species_tree.y,newick_string=newick_string.y)

#Remove duplicate model results
all_paml <- all_paml %>%
  unite(hog_treenum,model_num,col=hog_tree_model,remove=FALSE) %>%
  arrange(hog,treenum,model) %>%
  distinct(hog_tree_model,.keep_all=TRUE) %>%
  dplyr::select(-hog_tree_model)

#Read in HyPhy results
busted_res <- read_tsv("01_input_raw_paml_hyphy_results/busted_all_reruns_11.7.17.txt")
bsrel_res <- read.table("01_input_raw_paml_hyphy_results/bsrel_res_parsed_ratites_2017-11-01.txt",fill=TRUE) 
names(bsrel_res) <- c("class", "tree", "hog", "tsel.s", "nsel.s", "tsel.n", "nsel.n", "tnon", "nnon", "strict_branches", "nom_branches")

#Clean up bs-rel results
bsrel_res <- bsrel_res %>% tbl_df %>%
  mutate(total_sel.s = tsel.s + nsel.s, total_sel.n = tsel.n + nsel.n, totbranch = tnon + nnon + tsel.n + nsel.n, prop_sel.s = total_sel.s/totbranch, prop_sel.n = total_sel.n/totbranch) %>%
  dplyr::select(tree,hog,totbranch,total_sel.s,total_sel.n,prop_sel.s,prop_sel.n,strict_branches,nom_branches) %>%
  separate(tree,sep="ee",into=c("drop","treenum")) %>%
  unite(hog,treenum,col="hog_treenum") %>%
  dplyr::select(-drop) %>%
  filter(!is.na(total_sel.n))

#Clean up Busted results
busted_res <- busted_res %>%
  separate(treenum,sep = "ee",into=c("drop","treenum")) %>%
  unite(hog,treenum,col="hog_treenum") %>%
  dplyr::select(hog_treenum,pval_busted=pval,omega_busted=omega,weight_busted=weight) %>%
  filter(!is.na(omega_busted))

#Combine bs-rel and busted results:
hyphy_res <- busted_res %>%
  full_join(bsrel_res,by="hog_treenum")

#How many hogs?
length(unique(all_paml$hog))

#Which hogs are still missing a PAML model?
all_paml %>% group_by(hog_treenum) %>%
  summarise(tot_res = n()) %>%
  filter(tot_res != 6) %>%
  print(n=100)



#######################################################################################################################
#LRT and p-value calc for PAML results, and split into gene and species trees
#######################################################################################################################

#Function to conduct likelihood ratio tests between m1 and m2 models, and m7 and m8 models when given a hogid.
lrt_hog <- function(hogid,results){
  res_hog <- results %>% filter(hog==hogid)
  chisq_m1m2 <- chisq_m7m8 <- chisq_m2m2a <- chisq_m8m8a <- pval_m1m2 <- pval_m7m8 <- pval_m2m2a <- pval_m8m8a <- NA
  
  if ("1" %in% res_hog$model_num & "2" %in% res_hog$model_num){
    chisq_m1m2 <- 2*(res_hog[res_hog$model_num=="2","lnl"]-res_hog[res_hog$model_num=="1","lnl"]) %>% pull(lnl)
    pval_m1m2 <- 1-pchisq(chisq_m1m2,2)
  }
  
  if ("2" %in% res_hog$model_num & "2a" %in% res_hog$model_num){
    chisq_m2m2a <- 2*(res_hog[res_hog$model_num=="2","lnl"]-res_hog[res_hog$model_num=="2a","lnl"]) %>% pull(lnl)
    pval_m2m2a <- 1-pchisq(chisq_m2m2a,1)
  }
  
  if ("7" %in% res_hog$model_num & "8" %in% res_hog$model_num){
    chisq_m7m8 <- 2*2*(res_hog[res_hog$model_num=="8","lnl"]-res_hog[res_hog$model_num=="7","lnl"]) %>% pull(lnl)
    pval_m7m8 <- 1-pchisq(chisq_m7m8,2)
  }
  
  if ("8" %in% res_hog$model_num & "8a" %in% res_hog$model_num){
    chisq_m8m8a <- 2*2*(res_hog[res_hog$model_num=="8","lnl"]-res_hog[res_hog$model_num=="8a","lnl"]) %>% pull(lnl)
    pval_m8m8a <- 1-pchisq(chisq_m8m8a,1)
  }
  
  
  return(c(chisq_m1m2,pval_m1m2,chisq_m7m8,pval_m7m8,chisq_m2m2a,pval_m2m2a,chisq_m8m8a,pval_m8m8a))
}

#This function will take the results for a given hog, and extract the omega values and proportion of sites that fall under that omega for the m2 model.
omega_m2 <- function(hogid,results){
  res_hog <- results %>% filter(hog==hogid)
  prop_m2 <- omega_m2 <- NA
  
  if ("2" %in% res_hog$model_num){
    omega_m2_hog <- as.character(res_hog[res_hog$model_num=="2","omega"])
    omega_m2_hog <- strsplit(omega_m2_hog,split=",")[[1]]
    omega_m2_hog <- omega_m2_hog[length(omega_m2_hog)]
    omega_m2_hog <- strsplit(omega_m2_hog,split=":")[[1]]
    prop_m2 <- as.numeric(omega_m2_hog[1])
    omega_m2 <- as.numeric(omega_m2_hog[2])
  }
  
  return(c(prop_m2,omega_m2))
}

#This function will take the results for a given hog, and extract the omega values and proportion of sites that fall under that omega for the m8 model.
omega_m8 <- function(hogid,results){
  res_hog <- results %>% filter(hog==hogid)
  prop_m8 <- omega_m8 <- NA
  
  if ("8" %in% res_hog$model_num){
    omega_m8_hog <- as.character(res_hog[res_hog$model_num=="8","omega"])
    omega_m8_hog <- strsplit(omega_m8_hog,split=",")[[1]]
    omega_m8_hog <- omega_m8_hog[length(omega_m8_hog)]
    omega_m8_hog <- strsplit(omega_m8_hog,split=":")[[1]]
    prop_m8 <- as.numeric(omega_m8_hog[1])
    omega_m8 <- as.numeric(omega_m8_hog[2])
  }
  
  return(c(prop_m8,omega_m8))
}


#Split PAML results into gene and species trees, note that tree1 is the species tree - in cases where the gene tree and the species tree match (so there is only 1 tree, it is marked as "species_tree" is False.)
gene_trees <- paml_res_m0 %>%
  filter(species_tree == "False") %>%
  pull(hog_treenum)
# sp_trees <- paml_res_m0 %>%
#   separate(hog_treenum,into=c("hog","treenum")) %>%
#   filter(treenum == "1") %>%
#   unite(hog,treenum,col="hog_treenum") %>%
#   pull(hog_treenum)
sp_trees <- paml_res_m0 %>%
  filter(species_tree == "True") %>%
  pull(hog_treenum)

#Calculate p-values for all hogs between neutral and positive selection models (including FDR p-values)
#hogs_with_ids <- unique(paml_res_ncbi$hog)
all_hogs <- as.character(unique(all_paml$hog))
paml_pval_allgenes_gene <- matrix(nrow=length(all_hogs),ncol=13)
paml_pval_allgenes_sp <- matrix(nrow=length(all_hogs),ncol=13)

for (i in 1:length(all_hogs)){
  paml_pval_allgenes_gene[i,1:8] <- lrt_hog(all_hogs[i],all_paml[all_paml$hog_treenum %in% gene_trees,])
  paml_pval_allgenes_gene[i,9:10] <- try(omega_m2(all_hogs[i],all_paml[all_paml$hog_treenum %in% gene_trees,]))
  paml_pval_allgenes_gene[i,11:12] <- try(omega_m8(all_hogs[i],all_paml[all_paml$hog_treenum %in% gene_trees,]))
  if (length(grep(pattern=paste("^",paste(all_hogs[i],1,sep = "_"),sep=""),gene_trees,value=TRUE))>0){
    paml_pval_allgenes_gene[i,13] <- grep(pattern=paste("^",paste(all_hogs[i],1,sep = "_"),sep=""),gene_trees,value=TRUE)
  }
  if (length(grep(pattern=paste("^",paste(all_hogs[i],2,sep = "_"),sep=""),gene_trees,value=TRUE))>0){
    paml_pval_allgenes_gene[i,13] <- grep(pattern=paste("^",paste(all_hogs[i],2,sep = "_"),sep=""),gene_trees,value=TRUE)
  }
  
  paml_pval_allgenes_sp[i,1:8] <- lrt_hog(all_hogs[i],all_paml[all_paml$hog_treenum %in% sp_trees,])
  paml_pval_allgenes_sp[i,9:10] <- try(omega_m2(all_hogs[i],all_paml[all_paml$hog_treenum %in% sp_trees,]))
  paml_pval_allgenes_sp[i,11:12] <- try(omega_m8(all_hogs[i],all_paml[all_paml$hog_treenum %in% sp_trees,]))
  if (length(grep(pattern=paste("^",paste(all_hogs[i],1,sep = "_"),sep=""),sp_trees,value=TRUE))>0){
    paml_pval_allgenes_sp[i,13] <- grep(pattern=paste("^",paste(all_hogs[i],1,sep = "_"),sep=""),sp_trees,value=TRUE)
  }
  if (length(grep(pattern=paste("^",paste(all_hogs[i],2,sep = "_"),sep=""),sp_trees,value=TRUE))>0){
    paml_pval_allgenes_sp[i,13] <- grep(pattern=paste("^",paste(all_hogs[i],2,sep = "_"),sep=""),sp_trees,value=TRUE)
  }
  
}
paml_pval_allgenes_gene <- data.frame(paml_pval_allgenes_gene)
colnames(paml_pval_allgenes_gene) <- c("ChiSq_m1m2","PVal_m1m2","ChiSq_m7m8","PVal_m7m8","ChiSq_m2m2a","PVal_m2m2a","ChiSq_m8m8a","PVal_m8m8a","Prop_m2","Omega_m2","Prop_m8","Omega_m8","hog_treenum")
rownames(paml_pval_allgenes_gene) <- all_hogs

paml_pval_allgenes_sp <- data.frame(paml_pval_allgenes_sp)
colnames(paml_pval_allgenes_sp) <- c("ChiSq_m1m2","PVal_m1m2","ChiSq_m7m8","PVal_m7m8","ChiSq_m2m2a","PVal_m2m2a","ChiSq_m8m8a","PVal_m8m8a","Prop_m2","Omega_m2","Prop_m8","Omega_m8","hog_treenum")
rownames(paml_pval_allgenes_sp) <- all_hogs

paml_pval_allgenes_gene$FDRPval_m1m2 <- p.adjust(as.numeric(as.character(paml_pval_allgenes_gene$PVal_m1m2)),method="BH")
paml_pval_allgenes_gene$FDRPval_m2m2a <- p.adjust(as.numeric(as.character(paml_pval_allgenes_gene$PVal_m2m2a)),method="BH")
paml_pval_allgenes_gene$FDRPval_m7m8 <- p.adjust(as.numeric(as.character(paml_pval_allgenes_gene$PVal_m7m8)),method="BH")
paml_pval_allgenes_gene$FDRPval_m8m8a <- p.adjust(as.numeric(as.character(paml_pval_allgenes_gene$PVal_m8m8a)),method="BH")

paml_pval_allgenes_gene$hog <- rownames(paml_pval_allgenes_gene)

paml_pval_allgenes_sp$FDRPval_m1m2 <- p.adjust(as.numeric(as.character(paml_pval_allgenes_sp$PVal_m1m2)),method="BH")
paml_pval_allgenes_sp$FDRPval_m2m2a <- p.adjust(as.numeric(as.character(paml_pval_allgenes_sp$PVal_m2m2a)),method="BH")
paml_pval_allgenes_sp$FDRPval_m7m8 <- p.adjust(as.numeric(as.character(paml_pval_allgenes_sp$PVal_m7m8)),method="BH")
paml_pval_allgenes_sp$FDRPval_m8m8a <- p.adjust(as.numeric(as.character(paml_pval_allgenes_sp$PVal_m8m8a)),method="BH")

paml_pval_allgenes_sp$hog <- rownames(paml_pval_allgenes_sp)

#Need to turn all columns but last hog column from factors to numeric vectors
for (i in c(1:12,14:17)){
  paml_pval_allgenes_gene[,i] <- as.numeric(as.character(paml_pval_allgenes_gene[,i]))
  paml_pval_allgenes_sp[,i] <- as.numeric(as.character(paml_pval_allgenes_sp[,i]))
}

paml_pval_allgenes_gene <- paml_pval_allgenes_gene %>% tbl_df %>%
  mutate(hog_treenum = as.character(hog_treenum))
paml_pval_allgenes_sp <- paml_pval_allgenes_sp %>% tbl_df %>%
  mutate(hog_treenum = as.character(hog_treenum))

save(paml_pval_allgenes_gene,paml_pval_allgenes_sp,file="01_output_processed_data/paml_pvals_allgenes.RDat")

#######################################################################################################################
#FDR Hyphy results, combine datasets
#######################################################################################################################

#Calculate adusted p-values for Busted analyses
hyphy_res_gene <- hyphy_res %>%
  filter(hog_treenum %in% gene_trees) %>%
  mutate(FDRPval_busted = p.adjust(pval_busted,method="BH"))

hyphy_res_sp <- hyphy_res %>%
  filter(hog_treenum %in% sp_trees) %>%
  mutate(FDRPval_busted = p.adjust(pval_busted,method="BH"))

save(hyphy_res_gene,hyphy_res_sp,file="01_output_processed_data/hyphy_res_allgenes.RDat")

all_res_gene <- paml_pval_allgenes_gene %>% full_join(hyphy_res_gene,by="hog_treenum")
all_res_sp <- paml_pval_allgenes_sp %>% full_join(hyphy_res_sp,by="hog_treenum")

save(all_res_gene,all_res_sp,file="01_output_processed_data/all_res_allgenes.RDat")
