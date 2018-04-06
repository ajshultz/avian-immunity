setwd("~/Dropbox/BirdImmuneGeneEvolution/")

library(tidyverse)
library(RColorBrewer)
library(cowplot)

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

#Read in chicken to human ensembl ID conversion, join annotation dataset.
#ggal_hsap <- read_csv("05_input_mammal_data/galgal_ens_to_human.txt", col_names = c("ggalEnsId", "ensID", "confidence"))

#Join annotation dataset and ensembl ID conversion, create bip, pip and vip tibbles with only those genes that fall into each category.
#ggal_hsap <- full_join(mammal_annot, ggal_hsap, by="ensID") %>%
#  mutate(bip = !is.na(bip), vip=!is.na(vip), pip=!is.na(pip))
#ggal_annot <- ggal_hsap %>% dplyr::select(ggalEnsId, confidence, bip, pip, vip) %>% filter(!is.na(ggalEnsId)) %>% group_by(ggalEnsId) %>% summarize(bip = sum(bip)>0, vip=sum(vip)>0, pip=sum(pip)>0)

#Translate chicken ensembl numbers to NCBI numbers, join to annotation data and write to file
#ens_ncbi <- read_csv("05_input_mammal_data/galgal_ens_to_ncbi.txt", col_names = c("ggalEnsId", "entrezID", "name"))
#ggal_annot <- left_join(ggal_annot, ens_ncbi)
#write_csv(ggal_annot, "05_output_bird_mammal_comparison_results/ggal_alt_annotation.csv")

#Apply the vip, bip and pip annotations to the hogs in our dataset, then add them back into the full dataset
hog_annot <- all_res_sp_zf_hs %>%
  right_join(mammal_annot,by=c("hsapiens_homolog_ensembl_gene" = "ensID")) %>%
  group_by(hog) %>%
  select(hog,entrezgene,bip,vip,pip) %>%
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

#Add in data from primate specific testing, and combine both mammal datasets
primate<-read_tsv("05_input_mammal_data/primate-9sp-data.txt") %>%
  dplyr::select(ensID = Ensembl.Gene.ID, psr = PSR.total) %>%
  mutate(paml_sig = 1)

mammal_select <- full_join(enard, primate)

#mammal_select <- full_join(enard, primate) %>%
#  left_join(ggal_hsap)

#Remove missing data (replace by non-sig results, and filter non one-to-one mapping by choosing the minimum busted pvalue, etc.)
ggal_mammal_comp <- mammal_select %>%
  ungroup %>%
  dplyr::select(ggalEnsId, BUSTED:paml_sig) %>%
  filter(!is.na(ggalEnsId)) %>% 
  replace_na(list(BUSTED=1, bs_ct=0, psr=0, paml_sig=0)) %>% 
  group_by(ggalEnsId) %>% 
  summarize(bustedp = min(BUSTED, na.rm=T), bs_ct = max(bs_ct, na.rm=T), psr = max(psr, na.rm=T), paml_sig = max(paml_sig, na.rm=T))

#Add hog designations to mammal data, remove missing mammal data (replace by non-sig results, and filter non one-to-one mapping by choosing the minimum busted pvalue, etc.)
mammal_clean 
all_res_sp_zf_hs %>%
  select(ensembl_gene_id,hsapiens_homolog_ensembl_gene) %>%
  right_join(mammal_select,by=c("hsapiens_homolog_ensembl_gene" = "ensID")) %>%
  filter(!is.na(ensembl_gene_id)) %>%
  replace_na(list(BUSTED=1, bs_ct=0, psr=0, paml_sig=0)) %>% 
  group_by(ensembl_gene_id) %>% 
  summarize(bustedp = min(BUSTED, na.rm=T), bs_ct = max(bs_ct, na.rm=T), psr = max(psr, na.rm=T), paml_sig = max(paml_sig, na.rm=T))
  

#Add in ensembl IDs and write to file.
ggal_mammal_comp <- ggal_mammal_comp %>%
  left_join(ens_ncbi)

ggal_mammal_comp %>%
  write_csv("05_output_bird_mammal_comparison_results/ggal_mammal_comp.csv")

#Combine bird and mammal datasets, only keep BUSTED results to ensure direct comparisons with mammals.
birds<-all_res_sp_zf_hs %>%
  dplyr::select(entrezgene,hog,pval_busted:FDRPval_busted)

imm<-full_join(ggal_mammal_comp, ggal_annot) %>%
  filter(!is.na(entrezID)) %>%
  full_join(birds, by=c("entrezID" = "entrezgene")) %>% filter(!is.na(hog))

imm <- imm %>% distinct(hog, .keep_all=TRUE)

imm %>% filter(vip)

#############################################################################
###Use combined dataset to make inferences about similarity in selection###
#############################################################################

#Get -log10 pvalues.
imm<-imm %>% mutate(mammal_logp = -1 * log10(bustedp+2.22e-16), bird_logp = -1 * log10(pval_busted+2.22e-16))

#Plot correlations of p-values between birds and mammals, for each interaction-type.
vip_pval_plot <- imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>%
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


#imm %>% filter(mammal_logp > 8, bird_logp > 8) %>% dplyr::select(entrezID, vip) %>% print.data.frame

#boxplot(imm$bird_logp ~ imm$paml_sig)

#imm %>% dplyr::filter(paml_sig == 1, bird_logp < 4) %>% dplyr::select(name) %>% print.data.frame()

#boxplot(imm$bird_logp ~ imm$vip, notch=T)

#Create new q-values considering the genes in this overlap dataset
imm <- imm %>%
  mutate(mammal_q = p.adjust(bustedp, method="BH"), bird_q = p.adjust(pval_busted, method="BH"))





#Numbers of genes overlapping at q<0.1 - q<-0.0001 and fisher's exact tests for significance.
#imm <- imm %>% mutate(comp_busted = (mammal_q < 0.01) + (bird_q < 0.01))
#Create column to identify non-immune genes (not vips, bips, or pips)
imm <- imm %>% mutate(non_immune = as.logical(vip == FALSE & bip == FALSE & pip == FALSE))

qvals <- c(0.1,0.01,0.001,0.0001)

qval_res_list <- list()

for (i in 1:length(qvals)){
  
  comp_propsig <- matrix(nrow=5,ncol=12)
  comp_propsig[,1] <- c("non-immune","vips","bips","pips","all_genes")
  
  #Fisher's exact tests
  comp_propsig[1,2:9] <- imm %>% filter(!is.na(mammal_q), !is.na(bird_q)) %>% with(., table(non_immune, mammal_q < qvals[i] & bird_q < qvals[i])) %>% fisher.test %>% unlist
  comp_propsig[2,2:9] <- imm %>% filter(!is.na(mammal_q), !is.na(bird_q)) %>% with(., table(vip, mammal_q < qvals[i] & bird_q < qvals[i])) %>% fisher.test %>% unlist
  comp_propsig[3,2:9] <- imm %>% filter(!is.na(mammal_q), !is.na(bird_q)) %>% with(., table(bip, mammal_q < qvals[i] & bird_q < qvals[i])) %>% fisher.test %>% unlist
  comp_propsig[4,2:9] <- imm %>% filter(!is.na(mammal_q), !is.na(bird_q)) %>% with(., table(pip, mammal_q < qvals[i] & bird_q < qvals[i])) %>% fisher.test %>% unlist
  comp_propsig[5,2:9] <- NA
  
  #proportions selected in both
  comp_propsig[1,10:11] <- imm %>% filter(!is.na(mammal_q), !is.na(bird_q)) %>% filter(non_immune == TRUE) %>% with(., prop.table(table(mammal_q < qvals[i] & bird_q < qvals[i])))
  comp_propsig[2,10:11] <- imm %>% filter(!is.na(mammal_q), !is.na(bird_q)) %>% filter(vip == TRUE) %>% with(., prop.table(table(mammal_q < qvals[i] & bird_q < qvals[i])))
  comp_propsig[3,10:11] <- imm %>% filter(!is.na(mammal_q), !is.na(bird_q)) %>% filter(bip == TRUE) %>% with(., prop.table(table(mammal_q < qvals[i] & bird_q < qvals[i])))
  comp_propsig[4,10:11] <- imm %>% filter(!is.na(mammal_q), !is.na(bird_q)) %>% filter(pip == TRUE) %>% with(., prop.table(table(mammal_q < qvals[i] & bird_q < qvals[i])))
  comp_propsig[5,10:11] <- imm %>% filter(!is.na(mammal_q), !is.na(bird_q)) %>% with(., prop.table(table(mammal_q < qvals[i] & bird_q < qvals[i])))
  
  comp_propsig[,12] <- qvals[i]
  
  colnames(comp_propsig) <- c("class","p.value","conf.int1","conf.int2","estimated.odds.ratio","null.value.odds.ratio","alternative","method","data.name","prop_no","prop_sel_both","qval")
  
  #Clean up, select relevant columsn
  comp_propsig_clean <- comp_propsig %>%
    as.tibble %>%
    mutate(p.value=round(as.numeric(p.value),digits = 4),odds.ratio=round(as.numeric(estimated.odds.ratio),digits=2),prop_sel_both=round(as.numeric(prop_sel_both),digits=3), conf.lower=round(as.numeric(conf.int1),digits=3), conf.upper=round(as.numeric(conf.int2),digits=3)) %>%
    dplyr::select(class,qval,odds.ratio,conf.lower,conf.upper,p.value,prop_sel_both)

  
  qval_res_list[[i]] <- comp_propsig_clean
}

qval_res <- qval_res_list %>% bind_rows

write_csv(qval_res,path="05_output_bird_mammal_comparison_results/mammal_bird_prop_selected_test_allq_results.csv")


ggplot(qval_res,aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),prop_sel_both,fill=factor(class,levels=c("all_genes","non-immune","vips","bips","pips")))) + geom_bar(stat = "identity",position="dodge") +
  scale_x_discrete(labels=c("0.1","0.01","0.001","0.0001")) +
  scale_fill_discrete(name="class") +
  xlab("q-value") +
  ylab("proportion selected birds + mammals")
  


#Create barplot
ggplot(comp_propsig_clean, aes(factor(class,levels=c("all","vips","bips","pips")), prop_sel_both)) +
  geom_col(fill=c("gray30", "red2", "red2", "red2")) +
  theme_classic() +
  ylab("proportion selected") +
  xlab("class")
ggsave("05_output_bird_mammal_comparison_results/proportion_selected_birds_mammals.pdf")


###Compare proportions of genes selected in both mammals and birds at different q-value cutoffs
#Plot bird vs. mammal log pvals, color by sig category
#imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>% mutate(plotcol = paste0(mammal_q<0.01, bird_q<0.01)) %>% ggplot(aes(x=bird_logp, y=mammal_logp, col=plotcol)) + geom_point()


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
(imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>% with(., table(bird_q < 0.0001, mammal_q < 0.0001)) %>% chisq.test)$expected

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



sig_both <- imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>%
  filter(bird_q < 0.0001 & mammal_q < 0.0001) 








#Create column to ID hogs signficant in all site tests
hyphy_paml_res$sig <- rep("not_sig",nrow(hyphy_paml_res))
hyphy_paml_res[all_selected_hogs,"sig"] <- "sig"
