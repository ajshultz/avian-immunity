setwd("~/Dropbox/BirdImmuneGeneEvolution/")

library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(clusterProfiler)

#Load dataset with NCBI annotations
load("02_output_annotated_data/all_res_zf_hs.Rdat")

#Load BIP, PIP and VIP annotations
bip<-read_csv("05_input_mammal_data/bip_mammal.csv", col_names=c("ensID")) %>%
  mutate(bip=TRUE)
pip<-read_csv("05_input_mammal_data/pip_mammals.csv") %>% mutate(pip=TRUE) %>%
  dplyr::select(ensID=GeneID, pip)
vip<-read_csv("05_input_mammal_data/enard_vip.csv") %>% mutate(vip=TRUE) %>%
  dplyr::select(ensID=`Ensembl identifier`, vip)

#Combine datasets
mammal_annot <- bip %>%
  full_join(pip, by="ensID") %>%
  full_join(vip, by="ensID") %>%
  mutate(bip = !is.na(bip), vip=!is.na(vip), pip=!is.na(pip))

#Apply the vip, bip and pip annotations to the hogs in our dataset, then add them back into the full dataset
hog_annot <- all_res_sp_zf_hs %>%
  right_join(mammal_annot,by=c("ensembl_gene_id_hs" = "ensID")) %>%
  group_by(hog) %>%
  dplyr::select(hog,entrezgene,bip,vip,pip) %>%
  summarize(bip = sum(bip)>0, vip=sum(vip)>0, pip=sum(pip)>0)

all_res_sp_zf_hs <- all_res_sp_zf_hs %>%
  left_join(hog_annot,by="hog") %>%
  replace_na(list(vip = FALSE, bip = FALSE, pip = FALSE))

#Save hog vip, bip and pip annotations
write_csv(hog_annot, "05_output_bird_mammal_comparison_results/hog_alt_annotation.csv")


#Load mammal data
enard_orig<-read_csv("05_input_mammal_data/enard_hyphy.csv", col_names = T, guess_max=2000) %>% dplyr::rename(ensID = `Ensembl Gene ID`, BUSTED = `BUSTED P-value`)
enard <- enard_orig %>%
  gather(branch, propsel, Human:Elephant) %>%
  group_by(ensID, BUSTED) %>%
  summarize(bs_ct = sum(propsel>0)) %>%
  ungroup

#Create dataset with mammal data only annotated with immune gene categories:
mammals_only <- enard %>%
  left_join(mammal_annot) %>%
  replace_na(list(vip = FALSE, bip = FALSE, pip = FALSE)) %>%
  mutate(FDRPval_busted = p.adjust(BUSTED, method="BH"))

#Add in data from primate specific testing, and combine both mammal datasets
#primate<-read_tsv("05_input_mammal_data/primate-9sp-data.txt") %>%
#  dplyr::select(ensID = Ensembl.Gene.ID, psr = PSR.total) %>%
#  mutate(paml_sig = 1)

#mammal_select <- full_join(enard, primate)

#Add hog designations to mammal data, remove missing mammal data (replace by non-sig results, and filter non one-to-one mapping by choosing the minimum busted pvalue, etc.)
mammal_clean <- all_res_sp_zf_hs %>%
  ungroup() %>%
  dplyr::select(hog,ensembl_gene_id,ensembl_gene_id_hs) %>%
  right_join(enard,by=c("ensembl_gene_id_hs" = "ensID")) %>%
  filter(!is.na(hog)) %>%
  replace_na(list(BUSTED=1, bs_ct=0)) %>% 
  group_by(hog) %>% 
  summarize(bustedp = min(BUSTED, na.rm=T), bs_ct = max(bs_ct, na.rm=T))
  

#Combine bird and mammal datasets, only keep BUSTED results to ensure direct comparisons with mammals.
birds<-all_res_sp_zf_hs %>%
  mutate(sig_all = if_else(FDRPval_m1m2 < 0.05 & FDRPval_m2m2a < 0.05 & FDRPval_m7m8 < 0.05 & FDRPval_m8m8a < 0.05 & FDRPval_busted < 0.05, TRUE,FALSE)) %>%
  dplyr::select(entrezgene,entrezgene_hs,ensembl_gene_id_hs,hog,pval_busted:FDRPval_busted,sig_all,bip,vip,pip)

imm<-full_join(mammal_clean, birds) %>%
  filter(!is.na(hog)) %>%
  distinct(hog, .keep_all=TRUE) %>%
  filter(!is.na(bustedp), !is.na(pval_busted))


#Get -log10 pvalues and qvalues (only consider genes in both datasets)
imm <-imm %>%
  mutate(mammal_logp = -1 * log10(bustedp+2.22e-16), bird_logp = -1 * log10(pval_busted+2.22e-16)) %>%
  mutate(mammal_q = p.adjust(bustedp, method="BH"), bird_q = p.adjust(pval_busted, method="BH"))

write_csv(imm,path = "05_output_bird_mammal_comparison_results/bird_mammal_combined_dataset.csv")


###################### 
#What proportion of genes are under selection in both birds and mammals? Is there a significant overlap?
#####################

########
#Calculate numbers of genes overlapping at q<0.1 - q<-0.0001 and Fisher's exact tests for significance.
qvals <- c(0.1,0.01,0.001,0.0001)

comp_propsig <- matrix(nrow=4,ncol=11)

for (i in 1:length(qvals)){

  comp_propsig[i,1] <- qvals[i]
  
  #Fisher's exact tests
  comp_propsig[i,2:9] <- imm %>% with(., table(mammal_q < qvals[i], bird_q < qvals[i])) %>% fisher.test %>% unlist

  #number of genes in each category
  comp_propsig[i,10] <- imm %>% filter(!is.na(mammal_q), !is.na(bird_q)) %>% summarize(n()) %>% pull
  
  #number of genes significant in both in each category
  comp_propsig[i,11] <- imm %>% filter(!is.na(mammal_q), !is.na(bird_q)) %>%
    filter(mammal_q < qvals[i], bird_q < qvals[i]) %>%
    summarize(n()) %>% pull

  
  colnames(comp_propsig) <- c("qval","p.value","conf.int1","conf.int2","estimated.odds.ratio","null.value.odds.ratio","alternative","method","data.name","n.genes","n.sig.both")
 
}

#Clean up, select relevant columns
comp_propsig_clean <- comp_propsig %>%
  as.tibble %>%
  mutate(p.value=round(as.numeric(p.value),digits = 4),odds.ratio=round(as.numeric(estimated.odds.ratio),digits=2), conf.lower=round(as.numeric(conf.int1),digits=3), conf.upper=round(as.numeric(conf.int2),digits=3)) %>%
  dplyr::select(qval,n.genes,n.sig.both,odds.ratio,conf.lower,conf.upper,p.value)

write_csv(comp_propsig_clean,path="05_output_bird_mammal_comparison_results/mammal_bird_prop_selected_test_allq_overall_overlap.csv")

#Create vector of asterisks to include in plots:
comp_propsig_clean <- comp_propsig_clean %>%
  mutate(sig = case_when(
    is.na(p.value) ~ "",
    p.value > 0.05 ~ "",
    p.value <= 0.05 & p.value > 0.01 ~ "*",
    p.value <= 0.01 & p.value > 0.001 ~ "**",
    p.value <= 0.001 ~ "***"
  ))

comp_propsig_clean %>%
  ggplot(aes(factor(qval,levels=c("1e-04","0.001","0.01","0.1")),odds.ratio)) +
  geom_point(size=8,position=position_dodge(width=0.9),col="#44AA99") +
  geom_linerange(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),ymin=conf.lower,ymax=conf.upper),size=2,col="#44AA99") +
  geom_hline(aes(yintercept = 1),size=2,linetype="dashed",col="#882255") +
  scale_x_discrete(labels=c("0.0001","0.001","0.01","0.1")) +
  xlab("q-value") +
  ylab("odds ratio") +
  ylim(0,8) +
  theme(axis.title = element_text(size=24), axis.text = element_text(size=18))

ggsave(filename = "05_output_bird_mammal_comparison_results/mammal_bird_odds_ratio_selboth.pdf",width = 8, height=5)




#Read back in imm if running later
#imm <- read_csv("05_output_bird_mammal_comparison_results/bird_mammal_combined_dataset.csv")


imm <- imm %>%
  mutate(sig_birds_mammals =  case_when(
    bird_q<=qval & mammal_q>qvals[i] ~ "birds",
    bird_q>qval & mammal_q<=qvals[i] ~ "not_sig",
    bird_q<=qval & mammal_q<=qvals[i] ~ "birds_and_mammals",
    bird_q>qval & mammal_q>qvals[i] ~ "not_sig"),
    sig_mammals = if_else(mammal_q<=qvals[i],TRUE,FALSE),
    sig_birds = if_else(bird_q<=qvals[i],TRUE,FALSE))



###Which genes are under selection in both lineages?
#Looking at q = 0.01
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
  birds_mammals_df <- summary(birds_mammals_k)
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

#Read in previous bird results:
bird_sig_pathways <- read_csv("04_output_pathway_results/chicken_genetree_pathwayres_p1_q0.1.csv")

#Get sig bird pathways in bird mammal enrichment results
sig_pathways <- bird_sig_pathways %>%
  semi_join(enrich_res,by="Description") %>%
  pull(Description)

#Choose 10 pathways with the same nubmer of background sig genes, but not in these 10
not_sig_pathways <- c("Ubiquitin mediated proteolysis","Spliceosome","Focal adhesion","NOD-like receptor signaling pathway","Protein processing in endoplasmic reticulum","Regulation of actin cytoskeleton","Vascular smooth muscle contraction","Ribosome biogenesis in eukaryotes","Purine metabolism","Pyrimidine metabolism")

#pathway_colors <- c(brewer.pal(name = "PRGn",n=10),rep("grey",10))
pathway_colors<- c("#88CCEE","#999933","#882255","#332288","#DDCC77","#117733","#CC6677","#AA4499","#44AA99","#AA7744",rep("grey",10))
names(pathway_colors) <- c(sig_pathways,not_sig_pathways)


enrich_res %>%
  filter(Description %in% sig_pathways | Description %in% not_sig_pathways) %>%
  mutate(qvalue = round(qvalue,2)) %>%
  mutate(sig_qvals = if_else(qvalue<=0.2,1,0)) %>%
  ggplot(aes(log10(qval),enrichment,col=factor((Description),levels=c(sig_pathways,not_sig_pathways)))) +
  geom_line(size=2) +
  geom_point(size = 5) +
  ylab("fold enrichment") +
  scale_color_manual(values=pathway_colors,name="Pathway")

ggsave("05_output_bird_mammal_comparison_results/birds_mammals_enrichment_plot.pdf",width=9,height=7)

enrich_res %>%
  filter(Description %in% sig_pathways) %>%
  dplyr::select(qval,Description,enrichment,qvalue,BgRatio,GeneRatio)















#############################################################################
###Use combined dataset to make inferences about similarity in selection###
#############################################################################

#Plot correlations of p-values between birds and mammals, for each interaction-type.
vip_pval_plot <- imm %>%
  ggplot(aes(x=bird_logp, y=mammal_logp, col=vip, fill=vip)) +
  geom_point(size=0.5,alpha=0.5) + geom_smooth(method="lm") +
  theme_bw() +
  xlab("bird -log(p)") +
  ylab("mammal -log(p)") +
  theme(legend.position = "top")
bip_pval_plot <- imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>%
  ggplot(aes(x=bird_logp, y=mammal_logp, col=bip, fill=bip)) +
  geom_point(size=0.5,alpha=0.5) + geom_smooth(method="lm") +
  theme_bw() +
  xlab("bird -log(p)")+
  ylab("mammal -log(p)") +
  theme(legend.position = "top")
pip_pval_plot <- imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>%
  ggplot(aes(x=bird_logp, y=mammal_logp, col=pip, fill=pip)) +
  geom_point(size=0.5,alpha=0.5) + geom_smooth(method="lm") +
  theme_bw() +
  xlab("bird -log(p)")+
  ylab("mammal -log(p)") +
  theme(legend.position = "top")

plot_grid(vip_pval_plot, bip_pval_plot, pip_pval_plot,ncol=3)
ggsave(filename = "05_output_bird_mammal_comparison_results/mammal_bird_pval_plots.pdf",width = 8, height=3)


########
#Calculate numbers of genes in birds with q<0.1 - q<0.0001 and Fisher's exact testes for significance
birds_singles <- birds %>%
  distinct(hog,.keep_all=TRUE) %>%
  mutate(non_immune = as.logical(vip == FALSE & bip == FALSE & pip == FALSE)) %>%
  filter(hog %in% imm$hog)

qvals <- c(0.1,0.01,0.001,0.0001)
qval_res_list <- list()

for (i in 1:length(qvals)){
  
  comp_propsig <- matrix(nrow=5,ncol=13)
  comp_propsig[,1] <- c("non-immune","vips","bips","pips","all_genes")
  
  #Fisher's exact tests
  comp_propsig[1,2:9] <- birds_singles %>% with(., table(non_immune, FDRPval_busted < qvals[i])) %>% fisher.test %>% unlist
  comp_propsig[2,2:9] <- birds_singles %>% with(., table(vip, FDRPval_busted < qvals[i])) %>% fisher.test %>% unlist
  comp_propsig[3,2:9] <- birds_singles %>% with(., table(bip, FDRPval_busted < qvals[i])) %>% fisher.test %>% unlist
  comp_propsig[4,2:9] <- birds_singles %>% with(., table(pip, FDRPval_busted < qvals[i])) %>% fisher.test %>% unlist
  comp_propsig[5,2:9] <- NA
  
  #proportions selected in both
  comp_propsig[1,10:11] <- birds_singles %>% filter(non_immune == TRUE) %>% with(., prop.table(table(FDRPval_busted < qvals[i])))
  comp_propsig[2,10:11] <- birds_singles %>% filter(vip == TRUE) %>% with(., prop.table(table(FDRPval_busted < qvals[i])))
  comp_propsig[3,10:11] <- birds_singles %>% filter(bip == TRUE) %>% with(., prop.table(table(FDRPval_busted < qvals[i])))
  comp_propsig[4,10:11] <- birds_singles %>% filter(pip == TRUE) %>% with(., prop.table(table(FDRPval_busted < qvals[i])))
  comp_propsig[5,10:11] <- birds_singles %>% with(., prop.table(table(FDRPval_busted < qvals[i])))
  
  comp_propsig[,12] <- qvals[i]
  
  #number of genes in each category
  comp_propsig[1,13] <- birds_singles %>% filter(!is.na(FDRPval_busted)) %>% filter(non_immune == TRUE) %>% summarize(n()) %>% pull
  comp_propsig[2,13] <- birds_singles %>% filter(vip == TRUE) %>% summarize(n()) %>% pull
  comp_propsig[3,13] <- birds_singles %>% filter(bip == TRUE) %>% summarize(n()) %>% pull
  comp_propsig[4,13] <- birds_singles %>% filter(pip == TRUE) %>% summarize(n()) %>% pull
  comp_propsig[5,13] <- birds_singles %>% summarize(n()) %>% pull
  
  colnames(comp_propsig) <- c("class","p.value","conf.int1","conf.int2","estimated.odds.ratio","null.value.odds.ratio","alternative","method","data.name","prop_no","prop_sel_both","qval","n.genes")
  
  #Clean up, select relevant columns
  comp_propsig_clean <- comp_propsig %>%
    as.tibble %>%
    mutate(p.value=round(as.numeric(p.value),digits = 4),odds.ratio=round(as.numeric(estimated.odds.ratio),digits=2),prop_sel_both=round(as.numeric(prop_sel_both),digits=3), conf.lower=round(as.numeric(conf.int1),digits=3), conf.upper=round(as.numeric(conf.int2),digits=3)) %>%
    dplyr::select(class,qval,n.genes,odds.ratio,conf.lower,conf.upper,p.value,prop_sel_both)
  
  
  qval_res_list[[i]] <- comp_propsig_clean
}

qval_res_birds <- qval_res_list %>% bind_rows

write_csv(qval_res_birds,path="05_output_bird_mammal_comparison_results/birds_only_prop_selected_test_allq_results_joint_hogs_only.csv")


#######
#Calculate numbers of genes in mammals with q<0.1 - q<0.0001 and Fisher's exact testes for significance
mammals_only <- mammals_only %>%
  mutate(non_immune = as.logical(vip == FALSE & bip == FALSE & pip == FALSE)) %>%
#  filter(ensID %in% imm$ensembl_gene_id_hs)

qvals <- c(0.1,0.01,0.001,0.0001)
qval_res_list <- list()

for (i in 1:length(qvals)){
  
  comp_propsig <- matrix(nrow=5,ncol=13)
  comp_propsig[,1] <- c("non-immune","vips","bips","pips","all_genes")
  
  #Fisher's exact tests
  comp_propsig[1,2:9] <- mammals_only %>% with(., table(non_immune, FDRPval_busted < qvals[i])) %>% fisher.test %>% unlist
  comp_propsig[2,2:9] <- mammals_only %>% with(., table(vip, FDRPval_busted < qvals[i])) %>% fisher.test %>% unlist
  comp_propsig[3,2:9] <- mammals_only %>% with(., table(bip, FDRPval_busted < qvals[i])) %>% fisher.test %>% unlist
  comp_propsig[4,2:9] <- mammals_only %>% with(., table(pip, FDRPval_busted < qvals[i])) %>% fisher.test %>% unlist
  comp_propsig[5,2:9] <- NA
  
  #proportions selected in both
  comp_propsig[1,10:11] <- mammals_only %>% filter(non_immune == TRUE) %>% with(., prop.table(table(FDRPval_busted < qvals[i])))
  comp_propsig[2,10:11] <- mammals_only %>% filter(vip == TRUE) %>% with(., prop.table(table(FDRPval_busted < qvals[i])))
  comp_propsig[3,10:11] <- mammals_only %>% filter(bip == TRUE) %>% with(., prop.table(table(FDRPval_busted < qvals[i])))
  comp_propsig[4,10:11] <- mammals_only %>% filter(pip == TRUE) %>% with(., prop.table(table(FDRPval_busted < qvals[i])))
  comp_propsig[5,10:11] <- mammals_only %>% with(., prop.table(table(FDRPval_busted < qvals[i])))
  
  comp_propsig[,12] <- qvals[i]
  
  #number of genes in each category
  comp_propsig[1,13] <- mammals_only %>% filter(!is.na(FDRPval_busted)) %>% filter(non_immune == TRUE) %>% summarize(n()) %>% pull
  comp_propsig[2,13] <- mammals_only %>% filter(vip == TRUE) %>% summarize(n()) %>% pull
  comp_propsig[3,13] <- mammals_only %>% filter(bip == TRUE) %>% summarize(n()) %>% pull
  comp_propsig[4,13] <- mammals_only %>% filter(pip == TRUE) %>% summarize(n()) %>% pull
  comp_propsig[5,13] <- mammals_only %>% summarize(n()) %>% pull
  
  colnames(comp_propsig) <- c("class","p.value","conf.int1","conf.int2","estimated.odds.ratio","null.value.odds.ratio","alternative","method","data.name","prop_no","prop_sel_both","qval","n.genes")
  
  #Clean up, select relevant columns
  comp_propsig_clean <- comp_propsig %>%
    as.tibble %>%
    mutate(p.value=round(as.numeric(p.value),digits = 4),odds.ratio=round(as.numeric(estimated.odds.ratio),digits=2),prop_sel_both=round(as.numeric(prop_sel_both),digits=3), conf.lower=round(as.numeric(conf.int1),digits=3), conf.upper=round(as.numeric(conf.int2),digits=3)) %>%
    dplyr::select(class,qval,n.genes,odds.ratio,conf.lower,conf.upper,p.value,prop_sel_both)
  
  
  qval_res_list[[i]] <- comp_propsig_clean
}

qval_res_mammals <- qval_res_list %>% bind_rows

write_csv(qval_res_mammals,path="05_output_bird_mammal_comparison_results/mammals_only_prop_selected_test_allq_results_joint_hogs_only.csv")

#write_csv(qval_res_mammals,path="05_output_bird_mammal_comparison_results/mammals_only_prop_selected_test_allq_results.csv")


########
#Calculate numbers of genes overlapping at q<0.1 - q<-0.0001 and Fisher's exact tests for significance.
#Create column to identify non-immune genes (not vips, bips, or pips)
imm <- imm %>% mutate(non_immune = as.logical(vip == FALSE & bip == FALSE & pip == FALSE))

qvals <- c(0.1,0.01,0.001,0.0001)

qval_res_list <- list()

for (i in 1:length(qvals)){
  
  comp_propsig <- matrix(nrow=5,ncol=13)
  comp_propsig[,1] <- c("non-immune","vips","bips","pips","all_genes")
  
  #Fisher's exact tests
  comp_propsig[1,2:9] <- imm %>% with(., table(non_immune, mammal_q < qvals[i] & bird_q < qvals[i])) %>% fisher.test %>% unlist
  comp_propsig[2,2:9] <- imm %>% with(., table(vip, mammal_q < qvals[i] & bird_q < qvals[i])) %>% fisher.test %>% unlist
  comp_propsig[3,2:9] <- imm %>% with(., table(bip, mammal_q < qvals[i] & bird_q < qvals[i])) %>% fisher.test %>% unlist
  comp_propsig[4,2:9] <- imm %>% with(., table(pip, mammal_q < qvals[i] & bird_q < qvals[i])) %>% fisher.test %>% unlist
  comp_propsig[5,2:9] <- NA
  
  #proportions selected in both
  comp_propsig[1,10:11] <- imm %>% filter(non_immune == TRUE) %>% with(., prop.table(table(mammal_q < qvals[i] & bird_q < qvals[i])))
  comp_propsig[2,10:11] <- imm %>% filter(vip == TRUE) %>% with(., prop.table(table(mammal_q < qvals[i] & bird_q < qvals[i])))
  comp_propsig[3,10:11] <- imm %>% filter(bip == TRUE) %>% with(., prop.table(table(mammal_q < qvals[i] & bird_q < qvals[i])))
  comp_propsig[4,10:11] <- imm %>% filter(pip == TRUE) %>% with(., prop.table(table(mammal_q < qvals[i] & bird_q < qvals[i])))
  comp_propsig[5,10:11] <- imm %>% with(., prop.table(table(mammal_q < qvals[i] & bird_q < qvals[i])))
  
  comp_propsig[,12] <- qvals[i]
  
  #number of genes in each category
  comp_propsig[1,13] <- imm %>% filter(!is.na(mammal_q), !is.na(bird_q)) %>% filter(non_immune == TRUE) %>% summarize(n()) %>% pull
  comp_propsig[2,13] <- imm %>% filter(vip == TRUE) %>% summarize(n()) %>% pull
  comp_propsig[3,13] <- imm %>% filter(bip == TRUE) %>% summarize(n()) %>% pull
  comp_propsig[4,13] <- imm %>% filter(pip == TRUE) %>% summarize(n()) %>% pull
  comp_propsig[5,13] <- imm %>% summarize(n()) %>% pull
  
  colnames(comp_propsig) <- c("class","p.value","conf.int1","conf.int2","estimated.odds.ratio","null.value.odds.ratio","alternative","method","data.name","prop_no","prop_sel_both","qval","n.genes")
  
  #Clean up, select relevant columns
  comp_propsig_clean <- comp_propsig %>%
    as.tibble %>%
    mutate(p.value=round(as.numeric(p.value),digits = 4),odds.ratio=round(as.numeric(estimated.odds.ratio),digits=2),prop_sel_both=round(as.numeric(prop_sel_both),digits=3), conf.lower=round(as.numeric(conf.int1),digits=3), conf.upper=round(as.numeric(conf.int2),digits=3)) %>%
    dplyr::select(class,qval,n.genes,odds.ratio,conf.lower,conf.upper,p.value,prop_sel_both)

  
  qval_res_list[[i]] <- comp_propsig_clean
}

qval_res <- qval_res_list %>% bind_rows

write_csv(qval_res,path="05_output_bird_mammal_comparison_results/mammal_bird_prop_selected_test_allq_results.csv")

#Create vector of asterisks to include in plots:
qval_res <- qval_res %>%
  mutate(sig = case_when(
    is.na(p.value) ~ "",
    p.value > 0.05 ~ "",
    p.value <= 0.05 & p.value > 0.01 ~ "*",
    p.value <= 0.01 & p.value > 0.001 ~ "**",
    p.value <= 0.001 ~ "***"
  ))


prop_sel_plot <- qval_res %>%
  filter(class != "all_genes") %>%
  ggplot(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),prop_sel_both,fill=factor(class,levels=c("all_genes","non-immune","vips","bips","pips")))) +
  geom_bar(stat = "identity",position=position_dodge()) +
  geom_text(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),label=sig),position=position_dodge(width=0.9),size=10) +
  scale_x_discrete(labels=c("0.1","0.01","0.001","0.0001")) +
  scale_fill_discrete(name="class") +
  xlab("q-value") +
  ylab("proportion selected birds + mammals") +
  labs(subtitle="* = p < 0.05, ** = p < 0.01, *** = p < 0.001")
  
odds_ratio_plot <- qval_res %>%
  filter(class != "all_genes") %>%
  ggplot(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),odds.ratio,col=factor(class,levels=c("all_genes","non-immune","vips","bips","pips")))) +
  geom_point(size=4,position=position_dodge(width=0.9)) +
  geom_linerange(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),ymin=conf.lower,ymax=conf.upper,col=factor(class,levels=c("all_genes","non-immune","vips","bips","pips"))),position=position_dodge(width=0.9),size=1.5) +
  geom_hline(aes(yintercept = 1)) +
  scale_x_discrete(labels=c("0.1","0.01","0.001","0.0001")) +
  scale_color_discrete(name="class") +
  xlab("q-value") +
  ylab("odds ratio") +
  scale_y_continuous(labels=function(label) sprintf('%6.0f', label))

plot_grid(prop_sel_plot, odds_ratio_plot, ncol=1)
ggsave(filename = "05_output_bird_mammal_comparison_results/mammal_bird_prop_sel_both.pdf",width = 10, height=10)



#Logistic regression for mammal sig genes given bird sig genes, vip (etc), and the interaction between them

#q values to consider
qvals <- c(0.1,0.01,0.001,0.0001)

#List to capture results across q values
qval_res_int_list <- list()

to_test <- c("vip","bip","pip")

for (i in 1:length(qvals)){
  
  comp_propsig_int <- matrix(nrow=(length(to_test)*3),ncol=7)
  comp_propsig_int[,1] <- rep(to_test,each=3)
  
  start = 1
  
  for (j in 1:length(to_test)){
    #Create a vector of true and false for that evidence type for each hog, and remove duplicates
    imm_test <- imm %>%
      distinct(hog,.keep_all=TRUE) %>%
      mutate(mammal_sig = if_else(mammal_q<qvals[i],1,0),
             bird_sig = if_else(bird_q<qvals[i],1,0))
    
    int_form <- as.formula(paste0("mammal_sig ~ bird_sig*",to_test[j]))
    
    gl_res <- glm(int_form, family="binomial",data=imm_test)
    
    comp_propsig_int[start:(start+2),2] <- rownames(summary(gl_res)$coefficient[2:4,])
    
    comp_propsig_int[start:(start+2),3:6] <- summary(gl_res)$coefficient[2:4,]
    
    start = start+3
  }
  comp_propsig_int[,7] <- qvals[i]
  
  colnames(comp_propsig_int) <- c("class","parameter","estimate","std.error","z.value","p.value","qval")
  
  #Clean up, select relevant columns
  comp_propsig_int_clean <- comp_propsig_int %>%
    as.tibble %>%
    mutate(p.value=round(as.numeric(p.value),digits = 4),
           estimate=round(as.double(estimate),digits=2),
           std.error=round(as.numeric(std.error),digits=2),
           z.value=round(as.numeric(z.value),digits=2)) %>%
    dplyr::select(qval,class,parameter,estimate,std.error,z.value,p.value)
  
  qval_res_int_list[[i]] <- comp_propsig_int_clean
}

qval_res_int <- qval_res_int_list %>% bind_rows

write_csv(qval_res_int,path="05_output_bird_mammal_comparison_results/bird_mammals_comparison_logregression_table.csv")


#Produce plots for visualizing proportions of different classes
sig_birds_list <- list()
sig_mammals_list <- list()

for (i in 1:length(qvals)){
  
  #Add a column for sig genes birds, and sig genes mammals
  imm <- imm %>%
    mutate(sig_birds_mammals =  case_when(
      bird_q<=qval & mammal_q>qvals[i] ~ "birds_only",
      bird_q>qval & mammal_q<=qvals[i] ~ "mammals_only",
      bird_q<=qval & mammal_q<=qvals[i] ~ "birds_and_mammals",
      bird_q>qval & mammal_q>qvals[i] ~ "neither"),
      sig_mammals = if_else(mammal_q<=qvals[i],TRUE,FALSE),
      sig_birds = if_else(bird_q<=qvals[i],TRUE,FALSE))
  
  #Playing around with some plotting
  sig_birds_list[[i]] <- imm %>%
    ggplot(aes(vip,fill=sig_birds)) +
    geom_bar(position="fill") +
    facet_grid(~sig_mammals) +
    ggtitle(qvals[i])
  
  sig_mammals_list[[i]] <- imm %>%
    ggplot(aes(bip,fill=sig_mammals)) +
    geom_bar(position="fill") +
    facet_grid(~sig_birds)+
    ggtitle(qvals[i])
}

plot_grid(sig_birds_list[[1]],sig_birds_list[[2]],sig_birds_list[[3]],sig_birds_list[[4]],nrow=2,ncol=2)
ggsave(filename = "05_output_bird_mammal_comparison_results/Birds_sig_mammals_vips.png",height=6,width=12)
ggsave(filename = "05_output_bird_mammal_comparison_results/Birds_sig_mammals_pips.png",height=6,width=12)
ggsave(filename = "05_output_bird_mammal_comparison_results/Birds_sig_mammals_bips.png",height=6,width=12)
plot_grid(sig_mammals_list[[1]],sig_mammals_list[[2]],sig_mammals_list[[3]],sig_mammals_list[[4]],nrow=2,ncol=2)
ggsave(filename = "05_output_bird_mammal_comparison_results/Mammals_sig_birds_vips.png",height=6,width=12)
ggsave(filename = "05_output_bird_mammal_comparison_results/Mammals_sig_birds_pips.png",height=6,width=12)
ggsave(filename = "05_output_bird_mammal_comparison_results/Mammals_sig_birds_bips.png",height=6,width=12)



imm %>%
  summarize(n_vips = )
#Creat stacked bar plots

pips_colors <- c("white",rgb(27,158,119,100,maxColorValue = 225),rgb(27,158,119,225,maxColorValue = 225))
names(pips_colors) <- c("not significant","mammals","mammals and birds")

vips_colors <- c("white",rgb(117,112,179,100,maxColorValue = 225),rgb(117,112,179,225,maxColorValue = 225))
names(vips_colors) <- c("not significant","mammals","mammals and birds")

bips_colors <- c("white",rgb(217,95,2,100,maxColorValue = 225),rgb(217,95,2,225,maxColorValue = 225))
names(bips_colors) <- c("not significant","mammals","mammals and birds")

legend_colors <- c(rgb(102,102,102,225,maxColorValue = 225),rgb(102,102,102,100,maxColorValue = 225))
names(legend_colors) <- c("mammals and birds","mammals")

imm <- imm %>%
  mutate(sig_birds_mammals =  case_when(
    bird_q<=qval & mammal_q>qvals[i] ~ "not significant",
    bird_q>qval & mammal_q<=qvals[i] ~ "mammals",
    bird_q<=qval & mammal_q<=qvals[i] ~ "mammals and birds",
    bird_q>qval & mammal_q>qvals[i] ~ "not significant"),
    sig_mammals = if_else(mammal_q<=qvals[i],TRUE,FALSE),
    sig_birds = if_else(bird_q<=qvals[i],TRUE,FALSE))

mammals_pip <- imm %>%
  ggplot(aes(pip,fill=factor(sig_birds_mammals,levels=c("not significant","mammals and birds","mammals")))) +
  geom_bar(position="fill") +
  scale_fill_manual(values = pips_colors,name="significant",guide=FALSE) +
  ylab("proportion")
mammals_bip <- imm %>%
  ggplot(aes(bip,fill=factor(sig_birds_mammals,levels=c("not significant","mammals and birds","mammals")))) +
  geom_bar(position="fill") +
  scale_fill_manual(values = bips_colors,name="significant",guide=FALSE) +
  ylab("proportion")
mammals_vip <- imm %>%
  ggplot(aes(vip,fill=factor(sig_birds_mammals,levels=c("not significant","mammals and birds","mammals")))) +
  geom_bar(position="fill") +
  scale_fill_manual(values = vips_colors,name="significant",guide=FALSE) +
  ylab("proportion") +
  guides(fill="none")
mammals_legend <- imm %>%
  filter(sig_birds_mammals != "not significant") %>%
  ggplot(aes(vip,fill=factor(sig_birds_mammals,levels=c("mammals and birds","mammals")))) +
  geom_bar(position="fill") +
  scale_fill_manual(values = legend_colors,name="significant") +
  ylab("proportion") +
  ylim(1,2) +
  guides(fill=guide_legend(override.aes = list(fill=legend_colors,order=2))) +
  theme(line = element_blank(),
        axis.text = element_blank(),
        title = element_blank())

plot_grid(mammals_vip,mammals_bip,mammals_pip,mammals_legend,nrow=1)


#Create birds plots
pips_colors <- c("white",rgb(27,158,119,100,maxColorValue = 225),rgb(27,158,119,225,maxColorValue = 225))
names(pips_colors) <- c("not significant","birds","mammals and birds")

vips_colors <- c("white",rgb(117,112,179,100,maxColorValue = 225),rgb(117,112,179,225,maxColorValue = 225))
names(vips_colors) <- c("not significant","birds","mammals and birds")

bips_colors <- c("white",rgb(217,95,2,100,maxColorValue = 225),rgb(217,95,2,225,maxColorValue = 225))
names(bips_colors) <- c("not significant","birds","mammals and birds")

legend_colors <- c(rgb(102,102,102,225,maxColorValue = 225),rgb(102,102,102,100,maxColorValue = 225))
names(legend_colors) <- c("mammals and birds","birds")

imm <- imm %>%
  mutate(sig_birds_mammals =  case_when(
    bird_q<=qval & mammal_q>qvals[i] ~ "birds",
    bird_q>qval & mammal_q<=qvals[i] ~ "not significant",
    bird_q<=qval & mammal_q<=qvals[i] ~ "mammals and birds",
    bird_q>qval & mammal_q>qvals[i] ~ "not significant"),
    sig_mammals = if_else(mammal_q<=qvals[i],TRUE,FALSE),
    sig_birds = if_else(bird_q<=qvals[i],TRUE,FALSE))

birds_pip <- imm %>%
  ggplot(aes(pip,fill=factor(sig_birds_mammals,levels=c("not significant","mammals and birds","birds")))) +
  geom_bar(position="fill") +
  scale_fill_manual(values = pips_colors,name="significant",guide=FALSE) +
  ylab("proportion")
birds_bip <- imm %>%
  ggplot(aes(bip,fill=factor(sig_birds_mammals,levels=c("not significant","mammals and birds","birds")))) +
  geom_bar(position="fill") +
  scale_fill_manual(values = bips_colors,name="significant",guide=FALSE) +
  ylab("proportion")
birds_vip <- imm %>%
  ggplot(aes(vip,fill=factor(sig_birds_mammals,levels=c("not significant","mammals and birds","birds")))) +
  geom_bar(position="fill") +
  scale_fill_manual(values = vips_colors,name="significant",guide=FALSE) +
  ylab("proportion") +
  guides(fill="none")
birds_legend <- imm %>%
  filter(sig_birds_mammals != "not significant") %>%
  ggplot(aes(vip,fill=factor(sig_birds_mammals,levels=c("mammals and birds","birds")))) +
  geom_bar(position="fill") +
  scale_fill_manual(values = legend_colors,name="significant") +
  ylab("proportion") +
  ylim(1,2) +
  guides(fill=guide_legend(override.aes = list(fill=legend_colors,order=2))) +
  theme(line = element_blank(),
        axis.text = element_blank(),
        title = element_blank())

plot_grid(birds_vip,birds_bip,birds_pip,birds_legend,nrow=1)

plot_grid(mammals_vip,mammals_bip,mammals_pip,mammals_legend,birds_vip,birds_bip,birds_pip,birds_legend,nrow=2)
ggsave("05_output_bird_mammal_comparison_results/birds_mammals_vips_bips_pips_sig_proportions.pdf",height=10,width=16)











###Compare proportions of genes selected in both mammals and birds at different q-value cutoffs
#Plot bird vs. mammal log pvals, color by sig category

propsig_qvals <- matrix(nrow=4,ncol=11)
propsig_qvals[,1] <- c(0.1,0.01,0.001,0.0001)

#Fisher's exact tests
propsig_qvals[1,2:9] <- imm %>% with(., table(mammal_q < 0.1, bird_q < 0.1)) %>% fisher.test %>% unlist
propsig_qvals[2,2:9] <-imm %>% with(., table(mammal_q < 0.01, bird_q < 0.01)) %>% fisher.test %>% unlist
propsig_qvals[3,2:9] <- imm %>% with(., table(mammal_q < 0.001, bird_q < 0.001)) %>% fisher.test %>% unlist
propsig_qvals[4,2:9] <- imm %>% with(., table(mammal_q < 0.0001, bird_q < 0.0001)) %>% fisher.test %>% unlist
#P-values
propsig_qvals[1,10:11] <- imm %>% with(., prop.table(table(mammal_q < 0.1 & bird_q < 0.1)))
propsig_qvals[2,10:11] <- imm %>% with(., prop.table(table(mammal_q < 0.01 & bird_q < 0.01)))
propsig_qvals[3,10:11] <- imm %>% with(., prop.table(table(mammal_q < 0.001 & bird_q < 0.001)))
propsig_qvals[4,10:11] <- imm %>% with(., prop.table(table(mammal_q < 0.0001 & bird_q < 0.0001)))

colnames(propsig_qvals) <- c("class","p.value","conf.int1","conf.int2","estimated.odds.ratio","null.value.odds.ratio","alternative","method","data.name","prop_no","prop_sel_both")

#Clean up, select relevant columsn
propsig_qvals_clean <- propsig_qvals %>%
  as.tibble %>%
  mutate(p.value=round(as.numeric(p.value),digits = 4),odds.ratio=round(as.numeric(estimated.odds.ratio),digits=2),prop_sel_both=round(as.numeric(prop_sel_both),digits=3), conf.lower=round(as.numeric(conf.int1),digits=3), conf.upper=round(as.numeric(conf.int2),digits=3)) %>%
  dplyr::select(class,odds.ratio,conf.lower,conf.upper,p.value,prop_sel_both)

imm %>% filter(mammal_q<0.0001,bird_q<0.0001) %>% print(n=40)



(imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>% with(., table(bird_q < 0.1, mammal_q < 0.1)) %>% chisq.test)$observed
(imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>% with(., table(bird_q < 0.01, mammal_q < 0.01)) %>% fisher.test)
(imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>% with(., table(bird_q < 0.001, mammal_q < 0.001)) %>% chisq.test)$expected
(imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>% with(., table(bird_q < 0.0001, mammal_q < 0.0001)) %>% fisher.test)

imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>% with(., table(bird_q < 0.1, mammal_q < 0.1))
imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>% with(., table(bird_q < 0.01, mammal_q < 0.01))
imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>% with(., table(bird_q < 0.001, mammal_q < 0.001)) 
imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>% with(., table(bird_q < 0.0001, mammal_q < 0.0001)) 

ftestres <- data.frame(fdr=c("10%", "1%", "0.1%", "0.01%"), or=c(1.856, 2.107, 3.794,4.239), lc=c(1.602,1.7099,2.745,2.845), uc=c(2.152,2.6,5.279,6.374))
ggplot(ftestres, aes(x=fdr, y=or)) + ylim(0,6.5) +
  geom_errorbar(aes(ymin=lc, ymax=uc), width=0.1, col="purple", size=1) +
  geom_point(col="purple", size=5) + theme_classic(base_size = 20) + geom_hline(yintercept = 1, linetype=3, col="red", size=1)

imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp), !is.na(entrezID)) %>% select(entrezID) %>% write_tsv("background_set.txt", col_names = F)
imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp), !is.na(entrezID), bird_q < 0.0001, mammal_q < 0.0001) %>% select(entrezID) %>% write_tsv("select_set.txt", col_names = F)
imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp), !is.na(entrezID), bird_q < 0.0001 | mammal_q < 0.0001) %>% select(entrezID) %>% write_tsv("background2_set.txt", col_names = F)

imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp), !is.na(entrezID)) %>% select(entrezID) %>% write_tsv("background_set.txt", col_names = F)
imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp), !is.na(entrezID), bird_q < 0.01, mammal_q < 0.01) %>% select(entrezID) %>% write_tsv("select_set.txt", col_names = F)
imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp), !is.na(entrezID), bird_q < 0.01 | mammal_q < 0.01) %>% select(entrezID) %>% write_tsv("background2_set.txt", col_names = F)

imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp), !is.na(entrezID), mammal_q < 0.01, bird_q > .99) %>% select(entrezID) %>% write_tsv("mammal_specific.txt", col_names = F)
imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp), !is.na(entrezID), mammal_q > 0.99, bird_q < 0.01) %>% select(entrezID) %>% write_tsv("bird_specific.txt", col_names = F)

