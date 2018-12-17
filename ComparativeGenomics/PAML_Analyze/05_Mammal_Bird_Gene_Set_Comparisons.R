setwd("~/Dropbox/BirdImmuneGeneEvolution/")

library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(clusterProfiler)

#Load dataset with NCBI annotations
load("02_output_annotated_data/all_res_zf_hs.Rdat")

#Load mammal data
enard_orig<-read_csv("05_input_mammal_data/enard_hyphy.csv", col_names = T, guess_max=2000) %>% dplyr::rename(ensID = `Ensembl Gene ID`, BUSTED = `BUSTED P-value`)

#Clean up enard data
enard <- enard_orig %>%
  gather(branch, propsel, Human:Elephant) %>%
  group_by(ensID, BUSTED) %>%
  summarize(bs_ct = sum(propsel>0)) %>%
  ungroup %>%
  mutate(FDRPval_busted = p.adjust(BUSTED, method="BH"))


#Add hog designations to mammal data, remove missing mammal data (replace by non-sig results, and filter non one-to-one mapping by choosing the minimum busted pvalue, etc.)
mammal_clean <- all_res_sp_zf_hs %>%
  ungroup() %>%
  dplyr::select(hog,ensembl_gene_id,ensembl_gene_id_hs) %>%
  right_join(enard,by=c("ensembl_gene_id_hs" = "ensID")) %>%
  filter(!is.na(hog)) %>%
  replace_na(list(BUSTED=1, bs_ct=0)) %>% 
  group_by(hog) %>% 
  summarize(bustedp = min(BUSTED, na.rm=T), bs_ct = max(bs_ct, na.rm=T))
  

#Create bird dataset with simplified results (BUSTED only to ensure direct comparisons with mammals), plus note all genes sig with all tests
birds<-all_res_sp_zf_hs %>%
  mutate(sig_all = if_else(FDRPval_m1m2 < 0.05 & FDRPval_m2m2a < 0.05 & FDRPval_m7m8 < 0.05 & FDRPval_m8m8a < 0.05 & FDRPval_busted < 0.05, TRUE,FALSE)) %>%
  dplyr::select(entrezgene,entrezgene_hs,ensembl_gene_id_hs,hog,Omega_m0,pval_busted:FDRPval_busted,sig_all,length)

#Combine bird and mammal datasets, get -log10 pvalues and qvalues (only consider genes in both datasets)
imm<-full_join(mammal_clean, birds) %>%
  filter(!is.na(hog)) %>%
  distinct(hog, .keep_all=TRUE) %>%
  filter(!is.na(bustedp), !is.na(pval_busted)) %>%
  mutate(mammal_logp = -1 * log10(bustedp+2.22e-16), bird_logp = -1 * log10(pval_busted+2.22e-16)) %>%
  mutate(mammal_q = p.adjust(bustedp, method="BH"), bird_q = p.adjust(pval_busted, method="BH"))

#Save output to csv file
write_csv(imm,path = "05_output_bird_mammal_comparison_results/bird_mammal_combined_dataset.csv")




###################### 
#What proportion of genes are under selection in both birds and mammals? Is there a significant overlap?
#####################

#Exclude 20% of genes with the lowest m0 omega values (as calculated in birds)

imm_no20 <- imm %>%
  mutate(m0_rank = percent_rank(Omega_m0)) %>%
  filter(m0_rank > 0.2)


########
#Calculate numbers of genes overlapping at q<0.1 - q<-0.0001 and Fisher's exact tests for significance. Perform both with all genes and with 20% most constrained genes removed.
qvals <- c(0.1,0.01,0.001,0.0001)

comp_propsig <- matrix(nrow=4,ncol=11)
comp_propsig_no20 <- matrix(nrow=4,ncol=11)

log_reg_res_list <- list()
log_reg_res_list_no20 <- list()

for (i in 1:length(qvals)){

  comp_propsig[i,1] <- qvals[i]
  comp_propsig_no20[i,1] <- qvals[i]
  
  #Fisher's exact tests
  comp_propsig[i,2:9] <- imm %>% with(., table(mammal_q < qvals[i], bird_q < qvals[i])) %>% fisher.test %>% unlist
  comp_propsig_no20[i,2:9] <- imm_no20 %>% with(., table(mammal_q < qvals[i], bird_q < qvals[i])) %>% fisher.test %>% unlist

  #number of genes in each category
  comp_propsig[i,10] <- imm %>% filter(!is.na(mammal_q), !is.na(bird_q)) %>% summarize(n()) %>% pull
  comp_propsig_no20[i,10] <- imm_no20 %>% filter(!is.na(mammal_q), !is.na(bird_q)) %>% summarize(n()) %>% pull
  
  #number of genes significant in both in each category
  comp_propsig[i,11] <- imm %>% filter(!is.na(mammal_q), !is.na(bird_q)) %>%
    filter(mammal_q < qvals[i], bird_q < qvals[i]) %>%
    summarize(n()) %>% pull
  comp_propsig_no20[i,11] <- imm_no20 %>% filter(!is.na(mammal_q), !is.na(bird_q)) %>%
    filter(mammal_q < qvals[i], bird_q < qvals[i]) %>%
    summarize(n()) %>% pull
  
  #Logistic regression to include the alignment length to identify whether the overlap in bird and mammal results are driven by alignment length
  sink(paste0("05_output_bird_mammal_comparison_results/bird_mammal_sig_align_length_logistic_regression_sig_",qvals[i],".txt"))
  imm_test <- imm %>%
    mutate(mammal_sig = if_else(mammal_q<qvals[i],1,0),
           bird_sig = if_else(bird_q<qvals[i],1,0))
  
  gl_res <- glm("mammal_sig ~ bird_sig*length", family="binomial",data=imm_test)
  print(summary(gl_res))
  sink()
  
  sink(paste0("05_output_bird_mammal_comparison_results/bird_mammal_sig_align_length_logistic_regression_no20_sig_",qvals[i],".txt"))
  imm_test <- imm_no20 %>%
    mutate(mammal_sig = if_else(mammal_q<qvals[i],1,0),
           bird_sig = if_else(bird_q<qvals[i],1,0))
  
  gl_res <- glm("mammal_sig ~ bird_sig*length", family="binomial",data=imm_test)
  print(summary(gl_res)$coeff)
  sink()
  
  imm_test <- imm %>%
    mutate(mammal_sig = if_else(mammal_q<qvals[i],1,0),
           bird_sig = if_else(bird_q<qvals[i],1,0))
  
  gl_res <- glm("mammal_sig ~ bird_sig*length", family="binomial",data=imm_test)
  log_reg_res_list[[i]] <- summary(gl_res)$coeff %>%
    as.data.frame %>%
    rownames_to_column %>%
    as.tibble %>%
    mutate(qvalue = qvals[i])
  
  imm_test <- imm_no20 %>%
    mutate(mammal_sig = if_else(mammal_q<qvals[i],1,0),
           bird_sig = if_else(bird_q<qvals[i],1,0))
  
  gl_res <- glm("mammal_sig ~ bird_sig*length", family="binomial",data=imm_test)
  log_reg_res_list_no20[[i]] <- summary(gl_res)$coeff %>%
    as.data.frame %>%
    rownames_to_column %>%
    as.tibble %>%
    mutate(qvalue = qvals[i])

  
  colnames(comp_propsig) <- c("qval","p.value","conf.int1","conf.int2","estimated.odds.ratio","null.value.odds.ratio","alternative","method","data.name","n.genes","n.sig.both")
  colnames(comp_propsig_no20) <- c("qval","p.value","conf.int1","conf.int2","estimated.odds.ratio","null.value.odds.ratio","alternative","method","data.name","n.genes","n.sig.both")
 
}

#Clean up, select relevant columns
comp_propsig_clean <- comp_propsig %>%
  as.tibble %>%
  mutate(p.value=round(as.numeric(p.value),digits = 4),odds.ratio=round(as.numeric(estimated.odds.ratio),digits=2), conf.lower=round(as.numeric(conf.int1),digits=3), conf.upper=round(as.numeric(conf.int2),digits=3)) %>%
  dplyr::select(qval,n.genes,n.sig.both,odds.ratio,conf.lower,conf.upper,p.value)

write_csv(comp_propsig_clean,path="05_output_bird_mammal_comparison_results/mammal_bird_prop_selected_test_allq_overall_overlap.csv")

comp_propsig_clean_no20 <- comp_propsig_no20 %>%
  as.tibble %>%
  mutate(p.value=round(as.numeric(p.value),digits = 4),odds.ratio=round(as.numeric(estimated.odds.ratio),digits=2), conf.lower=round(as.numeric(conf.int1),digits=3), conf.upper=round(as.numeric(conf.int2),digits=3)) %>%
  dplyr::select(qval,n.genes,n.sig.both,odds.ratio,conf.lower,conf.upper,p.value)

write_csv(comp_propsig_clean_no20,path="05_output_bird_mammal_comparison_results/mammal_bird_prop_selected_test_allq_overall_overlap_excluding20perc.csv")

log_reg_res <- log_reg_res_list %>%
  bind_rows
write_csv(log_reg_res,path="05_output_bird_mammal_comparison_results/bird_mammal_sig_align_length_logistic_regression_all_res.csv")
log_reg_res_no20 <- log_reg_res_list_no20 %>%
  bind_rows
write_csv(log_reg_res_no20,path="05_output_bird_mammal_comparison_results/bird_mammal_sig_align_length_logistic_regression_excluding20perc_all_res.csv")


#Create vector of asterisks to include in plots:
comp_propsig_clean <- comp_propsig_clean %>%
  mutate(sig = case_when(
    is.na(p.value) ~ "",
    p.value > 0.05 ~ "",
    p.value <= 0.05 & p.value > 0.01 ~ "*",
    p.value <= 0.01 & p.value > 0.001 ~ "**",
    p.value <= 0.001 ~ "***"
  ))

#Plot odds ratio plot alone
comp_propsig_clean %>%
  ggplot(aes(factor(qval,levels=c("1e-04","0.001","0.01","0.1")),odds.ratio)) +
  geom_point(size=8,position=position_dodge(width=0.9),col="#44AA99") +
  geom_linerange(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),ymin=conf.lower,ymax=conf.upper),size=2,col="#44AA99") +
  geom_hline(aes(yintercept = 1),size=2,linetype="dashed",col="black") +
  scale_x_discrete(labels=c("0.0001","0.001","0.01","0.1")) +
  xlab("q-value") +
  ylab("odds ratio") +
  ylim(0,8) +
  theme(axis.title = element_text(size=24), axis.text = element_text(size=18))
ggsave(filename = "05_output_bird_mammal_comparison_results/mammal_bird_odds_ratio_selboth.pdf",width = 8, height=5)

comp_propsig_clean_no20 %>%
  ggplot(aes(factor(qval,levels=c("1e-04","0.001","0.01","0.1")),odds.ratio)) +
  geom_point(size=8,position=position_dodge(width=0.9),col="black") +
  geom_linerange(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),ymin=conf.lower,ymax=conf.upper),size=2,col="black") +
  geom_hline(aes(yintercept = 1),size=2,linetype="dashed",col="black") +
  scale_x_discrete(labels=c("0.0001","0.001","0.01","0.1")) +
  xlab("q-value") +
  ylab("odds ratio") +
  ylim(0,8) +
  theme(axis.title = element_text(size=24), axis.text = element_text(size=18))
ggsave(filename = "05_output_bird_mammal_comparison_results/mammal_bird_odds_ratio_selboth_no20.pdf",width = 8, height=5)

######################################################################################## 
###Is there any patheway enrichment for which genes are under selection in both lineages?
######################################################################################## 

#Look at with different q values
qvals <- c(0.1,0.01,0.001,0.0001)
enrich_res <- list()

for (i in 1:length(qvals)){
  imm <- imm %>%
    mutate(sig_birds_mammals =  case_when(
      bird_q<=qvals[i] & mammal_q>qvals[i] ~ "birds",
      bird_q>qvals[i] & mammal_q<=qvals[i] ~ "mammals",
      bird_q<=qvals[i] & mammal_q<=qvals[i] ~ "birds_and_mammals",
      bird_q>qvals[i] & mammal_q>qvals[i] ~ "not_sig"),
      sig_mammals = if_else(mammal_q<=qvals[i],TRUE,FALSE),
      sig_birds = if_else(bird_q<=qvals[i],TRUE,FALSE))
  
  sig_both_genes <- imm %>%
    filter(sig_birds_mammals == "birds_and_mammals") %>%
    pull(entrezgene) %>%
    as.character
  
  sig_birds_genes <- imm %>%
    filter(sig_birds) %>%
    pull(entrezgene) %>%
    as.character
  
  birds_mammals_k <- enrichKEGG(sig_both_genes,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=sig_birds_genes,keyType="ncbi-geneid")
  birds_mammals_df <- as.data.frame(birds_mammals_k)
  birds_mammals_df[,"qval"] <- qvals[i]
  enrich_res[[i]] <- as.tibble(birds_mammals_df)
}

enrich_res <- bind_rows(enrich_res)

enrich_res <- enrich_res %>%
  separate(GeneRatio,into=c("sig_genes_pathway","sig_genes"),remove = F) %>%
  separate(BgRatio, into=c("bg_sig_genes_pathway","bg_sig_genes"),remove = F) %>%
  mutate(GeneRatio = as.numeric(sig_genes_pathway)/as.numeric(sig_genes)) %>%
  mutate(enrichment = (as.numeric(sig_genes_pathway)/as.numeric(sig_genes))/(as.numeric(bg_sig_genes_pathway)/as.numeric(bg_sig_genes)))

enrich_res %>%
  write_csv("05_output_bird_mammal_comparison_results/bird_mammal_pathway_enrichment.csv")

#Read in previous bird-only results:
bird_sig_pathways <- read_csv("04_output_pathway_results/chicken_genetree_pathwayres_p1_q0.1.csv")

#Get sig bird pathways in bird mammal enrichment results
sig_pathways <- bird_sig_pathways %>%
  semi_join(enrich_res,by="Description") %>%
  pull(Description)

#Choose 10 pathways with the same nubmer of background sig genes, but not in these 10
not_sig_pathways <- c("Ubiquitin mediated proteolysis","Spliceosome","Focal adhesion","NOD-like receptor signaling pathway","Protein processing in endoplasmic reticulum","Regulation of actin cytoskeleton","Vascular smooth muscle contraction","Ribosome biogenesis in eukaryotes","Purine metabolism","Pyrimidine metabolism")

not_sig_pathways <- enrich_res %>%
  distinct(Description) %>%
  filter(!(Description %in% sig_pathways)) %>%
  pull

#pathway_colors <- c(brewer.pal(name = "PRGn",n=10),rep("grey",10))
pathway_colors<- c("#88CCEE","#999933","#332288","#882255","#DDCC77","#117733","#CC6677","#AA4499","#44AA99","#AA7744",rep("grey",length(not_sig_pathways)))
names(pathway_colors) <- c(sig_pathways,not_sig_pathways)

#Identify significant q values
enrich_res <- enrich_res %>%
  mutate(qvalue = round(qvalue,2)) %>%
  mutate(sig_qvals = if_else(qvalue<=0.2,1,0)) 

#Plot enrichment values for pathways sig in birds colored, all other pathways in grey. Plot for each q-value and connect with lines
enrich_res %>%
  filter(Description %in% sig_pathways) %>%
  ggplot(aes(log10(qval),enrichment,col=factor((Description),levels=c(sig_pathways,not_sig_pathways)))) +
  geom_line(data=enrich_res %>% filter(Description %in% not_sig_pathways), aes(log10(qval),enrichment,group=Description), col="grey", size=2,alpha=0.5) +
  geom_point(data=enrich_res %>% filter(Description %in% not_sig_pathways), aes(log10(qval),enrichment,group=Description), col="grey", size=4, alpha=0.5) +
  geom_point(data=enrich_res %>% filter(sig_qvals==1),col="black",size=8) +
  geom_line(size=2,alpha=0.5) +
  geom_point(size = 4,alpha=1) +
  ylab("enrichment") +
  scale_color_manual(values=pathway_colors,name="Pathway") +
  scale_x_continuous(labels=c("0.0001","0.001","0.01","0.1")) 
ggsave("05_output_bird_mammal_comparison_results/birds_mammals_enrichment_plot.pdf",width=9,height=7)

#Which pathways are enriched for selected genes in birds and mammals with q < 0.2?
enrich_res %>%
  filter(qvalue < 0.2)

#Perform permutation test to get estimated enrichement values for Influenza A and Herpes simplex pathways if random sets of genes are chosen from bird sig genes to enrich.
#For Influenza A

#q<0.0001
qval <- 0.0001
enrich_res_perm_0.0001_list <- list()

imm_0.0001 <- imm %>%
  mutate(sig_birds_mammals =  case_when(
    bird_q<=0.0001 & mammal_q>0.0001 ~ "birds",
    bird_q>0.0001 & mammal_q<=0.0001 ~ "mammals",
    bird_q<=0.0001 & mammal_q<=0.0001 ~ "birds_and_mammals",
    bird_q>0.0001 & mammal_q>0.0001 ~ "not_sig"),
    sig_mammals = if_else(mammal_q<=0.0001,TRUE,FALSE),
    sig_birds = if_else(bird_q<=0.0001,TRUE,FALSE))

sig_both_genes_0.0001 <- imm_0.0001 %>%
  filter(sig_birds_mammals == "birds_and_mammals") %>%
  pull(entrezgene) %>%
  as.character

sig_birds_genes_0.0001 <- imm_0.0001 %>%
  filter(sig_birds) %>%
  pull(entrezgene) %>%
  as.character

for (i in 1:1000){
  #Pull random sets of genes the same size as the actual number under selection in both lineages
  sig_both_genes_0.0001_rand <- sample(sig_birds_genes_0.0001,length(sig_both_genes_0.0001),replace = F)
  
  #Perform enrichment, save entire dataframe to list (as tibble)
  birds_mammals_k <- enrichKEGG(sig_both_genes_0.0001_rand,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=sig_birds_genes_0.0001,keyType="ncbi-geneid")
  birds_mammals_df <- as.data.frame(birds_mammals_k)
  birds_mammals_df[,"qval"] <- 0.0001
  birds_mammals_df[,"rep"] <- i
  enrich_res_perm_0.0001_list[[i]] <- as.tibble(birds_mammals_df)
}

enrich_res_perm_0.0001 <- bind_rows(enrich_res_perm_0.0001_list)

#Calculate enrichment:
enrich_res_perm_0.0001 <- enrich_res_perm_0.0001 %>%
  separate(GeneRatio,into=c("sig_genes_pathway","sig_genes"),remove = F) %>%
  separate(BgRatio, into=c("bg_sig_genes_pathway","bg_sig_genes"),remove = F) %>%
  mutate(GeneRatio = as.numeric(sig_genes_pathway)/as.numeric(sig_genes)) %>%
  mutate(enrichment = (as.numeric(sig_genes_pathway)/as.numeric(sig_genes))/(as.numeric(bg_sig_genes_pathway)/as.numeric(bg_sig_genes)))


#q<0.001
qval <- 0.001
enrich_res_perm_0.001_list <- list()

imm_0.001 <- imm %>%
  mutate(sig_birds_mammals =  case_when(
    bird_q<=0.001 & mammal_q>0.001 ~ "birds",
    bird_q>0.001 & mammal_q<=0.001 ~ "mammals",
    bird_q<=0.001 & mammal_q<=0.001 ~ "birds_and_mammals",
    bird_q>0.001 & mammal_q>0.001 ~ "not_sig"),
    sig_mammals = if_else(mammal_q<=0.001,TRUE,FALSE),
    sig_birds = if_else(bird_q<=0.001,TRUE,FALSE))

sig_both_genes_0.001 <- imm_0.001 %>%
  filter(sig_birds_mammals == "birds_and_mammals") %>%
  pull(entrezgene) %>%
  as.character

sig_birds_genes_0.001 <- imm_0.001 %>%
  filter(sig_birds) %>%
  pull(entrezgene) %>%
  as.character

for (i in 1:1000){
  #Pull random sets of genes the same size as the actual number under selection in both lineages
  sig_both_genes_0.001_rand <- sample(sig_birds_genes_0.001,length(sig_both_genes_0.001),replace = F)
  
  #Perform enrichment, save entire dataframe to list (as tibble)
  birds_mammals_k <- enrichKEGG(sig_both_genes_0.001_rand,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=sig_birds_genes_0.001,keyType="ncbi-geneid")
  birds_mammals_df <- as.data.frame(birds_mammals_k)
  birds_mammals_df[,"qval"] <- 0.001
  birds_mammals_df[,"rep"] <- i
  enrich_res_perm_0.001_list[[i]] <- as.tibble(birds_mammals_df)
}

enrich_res_perm_0.001 <- bind_rows(enrich_res_perm_0.001_list)

#Calculate enrichment:
enrich_res_perm_0.001 <- enrich_res_perm_0.001 %>%
  separate(GeneRatio,into=c("sig_genes_pathway","sig_genes"),remove = F) %>%
  separate(BgRatio, into=c("bg_sig_genes_pathway","bg_sig_genes"),remove = F) %>%
  mutate(GeneRatio = as.numeric(sig_genes_pathway)/as.numeric(sig_genes)) %>%
  mutate(enrichment = (as.numeric(sig_genes_pathway)/as.numeric(sig_genes))/(as.numeric(bg_sig_genes_pathway)/as.numeric(bg_sig_genes)))



#q<0.01
qval <- 0.01
enrich_res_perm_0.01_list <- list()

imm_0.01 <- imm %>%
  mutate(sig_birds_mammals =  case_when(
    bird_q<=0.01 & mammal_q>0.01 ~ "birds",
    bird_q>0.01 & mammal_q<=0.01 ~ "mammals",
    bird_q<=0.01 & mammal_q<=0.01 ~ "birds_and_mammals",
    bird_q>0.01 & mammal_q>0.01 ~ "not_sig"),
    sig_mammals = if_else(mammal_q<=0.01,TRUE,FALSE),
    sig_birds = if_else(bird_q<=0.01,TRUE,FALSE))

sig_both_genes_0.01 <- imm_0.01 %>%
  filter(sig_birds_mammals == "birds_and_mammals") %>%
  pull(entrezgene) %>%
  as.character

sig_birds_genes_0.01 <- imm_0.01 %>%
  filter(sig_birds) %>%
  pull(entrezgene) %>%
  as.character

for (i in 1:1000){
  #Pull random sets of genes the same size as the actual number under selection in both lineages
  sig_both_genes_0.01_rand <- sample(sig_birds_genes_0.01,length(sig_both_genes_0.01),replace = F)
  
  #Perform enrichment, save entire dataframe to list (as tibble)
  birds_mammals_k <- enrichKEGG(sig_both_genes_0.01_rand,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=sig_birds_genes_0.01,keyType="ncbi-geneid")
  birds_mammals_df <- as.data.frame(birds_mammals_k)
  birds_mammals_df[,"qval"] <- 0.01
  birds_mammals_df[,"rep"] <- i
  enrich_res_perm_0.01_list[[i]] <- as.tibble(birds_mammals_df)
}

enrich_res_perm_0.01 <- bind_rows(enrich_res_perm_0.01_list)

#Calculate enrichment:
enrich_res_perm_0.01 <- enrich_res_perm_0.01 %>%
  separate(GeneRatio,into=c("sig_genes_pathway","sig_genes"),remove = F) %>%
  separate(BgRatio, into=c("bg_sig_genes_pathway","bg_sig_genes"),remove = F) %>%
  mutate(GeneRatio = as.numeric(sig_genes_pathway)/as.numeric(sig_genes)) %>%
  mutate(enrichment = (as.numeric(sig_genes_pathway)/as.numeric(sig_genes))/(as.numeric(bg_sig_genes_pathway)/as.numeric(bg_sig_genes)))



#q<0.1
qval <- 0.1
enrich_res_perm_0.1_list <- list()

imm_0.1 <- imm %>%
  mutate(sig_birds_mammals =  case_when(
    bird_q<=0.1 & mammal_q>0.1 ~ "birds",
    bird_q>0.1 & mammal_q<=0.1 ~ "mammals",
    bird_q<=0.1 & mammal_q<=0.1 ~ "birds_and_mammals",
    bird_q>0.1 & mammal_q>0.1 ~ "not_sig"),
    sig_mammals = if_else(mammal_q<=0.1,TRUE,FALSE),
    sig_birds = if_else(bird_q<=0.1,TRUE,FALSE))

sig_both_genes_0.1 <- imm_0.1 %>%
  filter(sig_birds_mammals == "birds_and_mammals") %>%
  pull(entrezgene) %>%
  as.character

sig_birds_genes_0.1 <- imm_0.1 %>%
  filter(sig_birds) %>%
  pull(entrezgene) %>%
  as.character

for (i in 1:1000){
  #Pull random sets of genes the same size as the actual number under selection in both lineages
  sig_both_genes_0.1_rand <- sample(sig_birds_genes_0.1,length(sig_both_genes_0.1),replace = F)
  
  #Perform enrichment, save entire dataframe to list (as tibble)
  birds_mammals_k <- enrichKEGG(sig_both_genes_0.1_rand,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=sig_birds_genes_0.1,keyType="ncbi-geneid")
  birds_mammals_df <- as.data.frame(birds_mammals_k)
  birds_mammals_df[,"qval"] <- 0.1
  birds_mammals_df[,"rep"] <- i
  enrich_res_perm_0.1_list[[i]] <- as.tibble(birds_mammals_df)
}

enrich_res_perm_0.1 <- bind_rows(enrich_res_perm_0.1_list)

#Calculate enrichment:
enrich_res_perm_0.1 <- enrich_res_perm_0.1 %>%
  separate(GeneRatio,into=c("sig_genes_pathway","sig_genes"),remove = F) %>%
  separate(BgRatio, into=c("bg_sig_genes_pathway","bg_sig_genes"),remove = F) %>%
  mutate(GeneRatio = as.numeric(sig_genes_pathway)/as.numeric(sig_genes)) %>%
  mutate(enrichment = (as.numeric(sig_genes_pathway)/as.numeric(sig_genes))/(as.numeric(bg_sig_genes_pathway)/as.numeric(bg_sig_genes)))


#Combine all permutation tests for each qval
enrich_res_perm_all <- bind_rows(enrich_res_perm_0.0001,enrich_res_perm_0.001,enrich_res_perm_0.01,enrich_res_perm_0.1)

#Save enrichment results
save(enrich_res_perm_all,file="05_output_bird_mammal_comparison_results/bird_mammal_enrichment_permutation_tests.Rdat")

#Get p-value for each bird-significant pathway, compared to permutation test distribution.
pval_perm <- function(perm_enrich_data, obs_enrich_data, qval_test, pathway){
  #perm_enrich_data = tibble of enrichment results from permutations
  #obs_enrich_data = tibble of enrichment results from actual test
  #qval_test = qvalue to test
  #pathway = pathway to test
  
  test_dist <- perm_enrich_data %>% filter(Description == pathway,qval == qval_test) %>% pull(enrichment)
  obs_val <- obs_enrich_data %>% filter(Description == pathway, qval == qval_test) %>% pull(enrichment)
  
  #Perform exact calculation
  1-(sum(obs_val>test_dist)+1)/(length(test_dist)+1)
}



#Calculate p-values from permutations, and convert to q values within each qval tested
enrich_res_perm_pvals <- enrich_res %>%
  mutate(pval=purrr::map2_dbl(.x=Description, .y=qval,.f=function(x,y) pval_perm(obs_enrich_data=enrich_res,perm_enrich_data=enrich_res_perm_all,pathway=x, qval_test=y))) %>%
  group_by(qval) %>%
  mutate(qval_perm=p.adjust(pval,method = "BH"))


#Plot permutation results alone
enrich_res_perm_all %>%
  dplyr::filter(Description %in% sig_pathways) %>%
  ggplot(aes(enrichment,fill=factor(Description,levels=sig_pathways))) +
  geom_density(alpha=0.5) +
  scale_fill_manual(values=pathway_colors,name="") +
  scale_color_manual(values=pathway_colors,name="") +
  geom_vline(data=enrich_res%>%dplyr::filter(Description %in% sig_pathways) ,aes(xintercept=enrichment,col=factor(Description,levels=sig_pathways)), size=2) +
  facet_grid(factor(Description,levels=sig_pathways)~qval) 
ggsave("05_output_bird_mammal_comparison_results/bird_mammal_enrichment_permutation_test_plot.pdf",width=10,height=8)


#Plot all together
#Plot odds ratio for use with other polots
odds_ratio_plot <- comp_propsig_clean %>%
  ggplot(aes(factor(qval,levels=c("1e-04","0.001","0.01","0.1")),odds.ratio)) +
  geom_point(size=4,position=position_dodge(width=0.9),col="black") +
  geom_linerange(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),ymin=conf.lower,ymax=conf.upper),size=2,col="black") +
  geom_hline(aes(yintercept = 1),size=2,linetype="dashed",col="black") +
  scale_x_discrete(labels=c("0.0001","0.001","0.01","0.1")) +
  xlab("q-value") +
  ylab("odds ratio") +
  ylim(0,7.5)

#Plot enrich res for use with other plots
all_enrich_plot <- enrich_res %>%
  filter(Description %in% sig_pathways) %>%
  ggplot(aes(log10(qval),enrichment,col=factor((Description),levels=c(sig_pathways,not_sig_pathways)))) +
  geom_line(data=enrich_res %>% filter(Description %in% not_sig_pathways), aes(log10(qval),enrichment,group=Description), col="grey", size=2,alpha=0.5) +
  geom_point(data=enrich_res %>% filter(Description %in% not_sig_pathways), aes(log10(qval),enrichment,group=Description), col="grey", size=4, alpha=0.5) +
  geom_point(data=enrich_res %>% filter(qvalue<.2),col="grey46",size=8) +
  geom_point(data=enrich_res %>% filter(qvalue<.1),col="black",size=8) +
  geom_line(size=2,alpha=0.5) +
  geom_point(size = 4,alpha=1) +
  ylab("enrichment") +
  scale_color_manual(values=pathway_colors,name="Pathway",guide=F) +
  scale_x_continuous(labels=c("0.0001","0.001","0.01","0.1"),name="q-value") +
  coord_cartesian(ylim=c(0,7.5))

perm_plot <- enrich_res_perm_all %>%
  dplyr::filter(Description %in% sig_pathways) %>%
  ggplot(aes(enrichment,fill=factor(Description,levels=sig_pathways))) +
  geom_density(alpha=0.5) +
  scale_fill_manual(values=pathway_colors,name="",guide=F) +
  scale_color_manual(values=pathway_colors,name="",guide=F) +
  #geom_vline(data=enrich_res_perm_pvals %>% filter(Description %in% sig_pathways,pval<0.05), aes(xintercept=enrichment,group=factor(Description,levels=sig_pathways)),size=3,color="grey") +
  #geom_vline(data=enrich_res_perm_pvals %>% filter(Description %in% sig_pathways,pval<0.01), aes(xintercept=enrichment,group=factor(Description,levels=sig_pathways)),size=3,color="black") +
  geom_vline(data=enrich_res_perm_pvals %>% filter(Description %in% sig_pathways,qval_perm<0.2), aes(xintercept=enrichment,group=factor(Description,levels=sig_pathways)),size=4,color="grey46") +
  geom_vline(data=enrich_res_perm_pvals %>% filter(Description %in% sig_pathways,qval_perm<0.1), aes(xintercept=enrichment,group=factor(Description,levels=sig_pathways)),size=4,color="black") +
  geom_vline(data=enrich_res%>%dplyr::filter(Description %in% sig_pathways) ,aes(xintercept=enrichment,col=factor(Description,levels=sig_pathways)), size=2) +
  scale_y_continuous(breaks=c(0,0.5,1)) +
  coord_cartesian(xlim=c(0,7)) +
  facet_grid(factor(Description,levels=sig_pathways)~qval) +
  theme(strip.text.y=element_blank(),
        axis.text.y=element_text(size=8))

pathway_legend <- enrich_res_perm_all %>%
  dplyr::filter(Description %in% sig_pathways) %>%
  ggplot(aes(enrichment,fill=factor(Description,levels=sig_pathways))) +
  geom_density(alpha=0.5) +
  scale_fill_manual(values=pathway_colors,name="") +
  scale_color_manual(values=pathway_colors,name="") +
  geom_vline(data=enrich_res%>%dplyr::filter(Description %in% sig_pathways) ,aes(xintercept=enrichment,col=factor(Description,levels=sig_pathways)), size=2) +
  coord_cartesian(xlim=c(8,8.5),ylim=1.5,1.8)+
  theme(line = element_blank(),
        axis.text = element_blank(),
        title = element_blank())

ggdraw() +
  draw_plot(odds_ratio_plot, x=0,y=0.6,width=0.5,height=0.4) +
  draw_plot(all_enrich_plot, x=0.5,y=0.6, width=0.5,height=0.4) +
  draw_plot(perm_plot,x=0,y=0, width=0.7,height=0.6) +
  draw_plot(pathway_legend,x=0.7,y=0,width=0.3,height=0.6) +
  draw_plot_label(label = c("A","B","C"), size=15,x=c(0,0.5,0),y=c(1,1,0.6))
ggsave("05_output_bird_mammal_comparison_results/bird_mammal_enrichment_all_plots_figure.pdf",width=13,height=8)
