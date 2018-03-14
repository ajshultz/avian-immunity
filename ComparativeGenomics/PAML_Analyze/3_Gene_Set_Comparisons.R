setwd("~/Dropbox/BirdImmuneGeneEvolution/")

#library(tidyverse)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)


#Load dataset created with 1_PAML_Hyphy_res_DataPrep.R
load("~/Dropbox/BirdImmuneGeneEvolution/PAML/processed_data/all_res_ncbi.Rdat")

#How many genes with all tests for gene trees:
all_res_gene_ncbi %>%
  filter(!is.na(pval_busted) & !is.na(PVal_m1m2) & !is.na(PVal_m2m2a) & !is.na(PVal_m7m8) & !is.na(PVal_m8m8a) & !is.na(total_sel.n)) %>%
  summarise(n())

#How many genes with all tests for species trees:
all_res_sp_ncbi %>%
  filter(!is.na(pval_busted) & !is.na(PVal_m1m2) & !is.na(PVal_m2m2a) & !is.na(PVal_m7m8) & !is.na(PVal_m8m8a) & !is.na(total_sel.n)) %>%
  summarise(n())

#Load BIP, PIP and VIP annotations
bip<-read_csv("Mammal_Comparisons/bip_mammal.csv", col_names=c("ensID")) %>%
  mutate(bip=TRUE)
pip<-read_csv("Mammal_Comparisons/pip_mammals.csv") %>% mutate(pip=TRUE) %>%
  dplyr::select(ensID=GeneID, pip)
vip<-read_csv("Mammal_Comparisons/enard_vip.csv") %>% mutate(vip=TRUE) %>%
  dplyr::select(ensID=`Ensembl identifier`, vip)

#Combine datasets
mammal_annot <- bip %>%
  full_join(pip, by="ensID") %>%
  full_join(vip, by="ensID")

#Read in chicken to human ensembl ID conversion, join annotation dataset.
ggal_hsap <- read_csv("Mammal_Comparisons/galgal_ens_to_human.txt", col_names = c("ggalEnsId", "ensID", "confidence"))

#Join annotation dataset and ensembl ID conversion, create bip, pip and vip tibbles with only those genes that fall into each category.
ggal_hsap <- full_join(mammal_annot, ggal_hsap, by="ensID") %>%
  mutate(bip = !is.na(bip), vip=!is.na(vip), pip=!is.na(pip))
ggal_annot <- ggal_hsap %>% dplyr::select(ggalEnsId, confidence, bip, pip, vip) %>% filter(!is.na(ggalEnsId)) %>% group_by(ggalEnsId) %>% summarize(bip = sum(bip)>0, vip=sum(vip)>0, pip=sum(pip)>0)

#Translate chicken ensembl numbers to NCBI numbers, join to annotation data and write to file
ens_ncbi <- read_csv("Mammal_Comparisons/galgal_ens_to_ncbi.txt", col_names = c("ggalEnsId", "entrezID", "name"))
ggal_annot <- left_join(ggal_annot, ens_ncbi)
write_csv(ggal_annot, "Mammal_Comparisons/ggal_alt_annotation.csv")

#Load mammal data
enard_orig<-read_csv("Mammal_Comparisons/enard_hyphy.csv", col_names = T, guess_max=2000) %>% dplyr::rename(ensID = `Ensembl Gene ID`, BUSTED = `BUSTED P-value`)
enard <- enard_orig %>%
  gather(branch, propsel, Human:Elephant) %>%
  group_by(ensID, BUSTED) %>%
  summarize(bs_ct = sum(propsel>0)) %>%
  ungroup

#Add in data from primate specific testing
primate<-read_tsv("Mammal_Comparisons/primate-9sp-data.txt") %>%
  dplyr::select(ensID = Ensembl.Gene.ID, psr = PSR.total) %>%
  mutate(paml_sig = 1)

mammal_select = full_join(enard, primate) %>%
  left_join(ggal_hsap)

#Remove missing data (replace by non-sig results, and filter non one-to-one mapping by choosing the minimum busted pvalue, etc.)
ggal_mammal_comp = mammal_select %>%
  ungroup %>%
  dplyr::select(ggalEnsId, BUSTED:paml_sig) %>%
  filter(!is.na(ggalEnsId)) %>% 
  replace_na(list(BUSTED=1, bs_ct=0, psr=0, paml_sig=0)) %>% 
  group_by(ggalEnsId) %>% 
  summarize(bustedp = min(BUSTED, na.rm=T), bs_ct = max(bs_ct, na.rm=T), psr = max(psr, na.rm=T), paml_sig = max(paml_sig, na.rm=T))

#Add in ensembl IDs and write to file.
ggal_mammal_comp <- ggal_mammal_comp %>%
  left_join(ens_ncbi)

ggal_mammal_comp %>%
  write_csv("Mammal_Comparisons/ggal_mammal_comp.csv")

#Combine bird and mammal datasets, only keep BUSTED results to ensure direct cmoparisons with mammals.
birds<-all_res_sp_ncbi %>%
  dplyr::select(entrezgene,hog,pval_busted:FDRPval_busted)

imm<-full_join(ggal_mammal_comp, ggal_annot) %>%
  filter(!is.na(entrezID)) %>%
  full_join(., birds, by=c("entrezID" = "entrezgene")) %>% filter(!is.na(hog))

#lets see if there is overlap in sig

#Get -log10 pvalues.
imm<-imm %>% mutate(mammal_logp = -1 * log10(bustedp+2.22e-16), bird_logp = -1 * log10(pval_busted+2.22e-16))

#Plot correlations of p-values between birds and mammals, for each interaction-type.
imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>%
  ggplot(aes(x=bird_logp, y=mammal_logp, col=vip)) + geom_point() + geom_smooth(method="lm")
imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>%
  ggplot(aes(x=bird_logp, y=mammal_logp, col=bip)) + geom_point() + geom_smooth(method="lm")
imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>%
  ggplot(aes(x=bird_logp, y=mammal_logp, col=pip)) + geom_point() + geom_smooth(method="lm")

#imm %>% filter(mammal_logp > 8, bird_logp > 8) %>% dplyr::select(entrezID, vip) %>% print.data.frame

#boxplot(imm$bird_logp ~ imm$paml_sig)

#imm %>% dplyr::filter(paml_sig == 1, bird_logp < 4) %>% dplyr::select(name) %>% print.data.frame()

#boxplot(imm$bird_logp ~ imm$vip, notch=T)

#Create new q-values considering the genes in this overlap dataset
imm <- imm %>% mutate(mammal_q = p.adjust(bustedp, method="BH"), bird_q = p.adjust(pval_busted, method="BH"))
table(imm$mammal_q < 0.01, imm$bird_q < 0.01) %>% fisher.test

#Numbers of genes overlapping at different q-values and fisher's exact tests for significance.
imm <- imm %>% mutate(comp_busted = (mammal_q < 0.01) + (bird_q < 0.01))

imm %>% filter(vip == FALSE, pip == FALSE) %>% with(., table(mammal_q < 0.001, bird_q < 0.001))
imm %>% filter(vip == TRUE | pip == TRUE) %>% with(., table(mammal_q < 0.001, bird_q < 0.001))

imm %>% with(., table(vip, mammal_q < 0.01 & bird_q < 0.01)) %>% fisher.test
imm %>% with(., table(bip, mammal_q < 0.01 & bird_q < 0.01)) %>% fisher.test
imm %>% with(., table(pip, mammal_q < 0.01 & bird_q < 0.01)) %>% fisher.test

imm %>% filter(vip == FALSE, pip == FALSE, bip == FALSE) %>% with(., prop.table(table(mammal_q < 0.01 & bird_q < 0.01)))
imm %>% filter(vip == TRUE) %>% with(., prop.table(table(mammal_q < 0.01 & bird_q < 0.01)))
imm %>% filter(pip == TRUE) %>% with(., prop.table(table(mammal_q < 0.01 & bird_q < 0.01)))
imm %>% filter(bip == TRUE) %>% with(., prop.table(table(mammal_q < 0.01 & bird_q < 0.01)))

#Need to create this dataframe with above results
prop_overlap_data<-data.frame(class=factor(c("non-immune", "bip", "vip", "pip"), levels=c("non-immune", "bip", "vip", "pip")), prop=c(0.02710997, 0.04973822, 0.04964539, 0.08715596))
ggplot(prop_overlap_data, aes(class, prop)) + geom_col(fill=c("gray30", "red2", "red2", "red2")) + theme_classic()

#genome graph
library(RColorBrewer)
genomes<-data.frame(group = c("mammals", "birds", "fishes", "amphibians", "reptiles"), count=c(24+20+20+45, 61, 46, 13, 3))
pie(genomes$count, labels=genomes$group, col=brewer.pal(5, "Set1"), cex=2)

#add a color for plot

imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>% mutate(plotcol = paste0(mammal_q<0.01, bird_q<0.01)) %>% ggplot(aes(x=bird_logp, y=mammal_logp, col=plotcol)) + geom_point()

imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>% with(., table(mammal_q < 0.01))
imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>% with(., table(bird_q < 0.01))

(imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>% with(., table(bird_q < 0.1, mammal_q < 0.1)) %>% chisq.test)$observed
(imm %>% filter(!is.na(mammal_logp), !is.na(bird_logp)) %>% with(., table(bird_q < 0.01, mammal_q < 0.01)) %>% chisq.test)$expected
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




#How many hogs are selected according to BUSTED and the PAML (m2 vs m2a and m8 vs m8a)
all_selected_hogs <- hyphy_paml_res[hyphy_paml_res[,"FDRPval_m2m2a"]<0.05 & hyphy_paml_res[,"FDRPval_m8m8a"]<0.05 & hyphy_paml_res[,"FDRPval_m1m2"]<0.05 & hyphy_paml_res[,"FDRPval_m7m8"]<0.05 & hyphy_paml_res[,"FDRPval_busted"]<0.05,"hog"]
all_selected_hogs <- as.character(all_selected_hogs[!is.na(all_selected_hogs)])
length(all_selected_hogs)

#How many branches are selected by BS-Rel for HOGs selected by all tests
hist(hyphy_paml_res[all_selected_hogs,"total_sel.s"],breaks=20)
hist(hyphy_paml_res[all_selected_hogs,"total_sel.n"],breaks=20)
summary(hyphy_paml_res[all_selected_hogs,"total_sel.s"])
summary(hyphy_paml_res[all_selected_hogs,"total_sel.n"])



#Create column to ID hogs signficant in all site tests
hyphy_paml_res$sig <- rep("not_sig",nrow(hyphy_paml_res))
hyphy_paml_res[all_selected_hogs,"sig"] <- "sig"
