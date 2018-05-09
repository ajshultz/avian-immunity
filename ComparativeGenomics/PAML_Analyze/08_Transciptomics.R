setwd("~/Dropbox/BirdImmuneGeneEvolution")
library(tidyverse)
library(biomaRt)
library(cowplot)

#Read in bird and mammal selection data
imm <- read_csv("05_output_bird_mammal_comparison_results/bird_mammal_combined_dataset.csv")

#Read in list of files
sleuth_files <- list.files("08_inputs_transcriptomics/sleuth/")

#Grab bioproject IDs from "prep" file
bioproj <- grep(sleuth_files,pattern=".prep",value = TRUE) %>%
  purrr::map_chr(.x=.,.f =function(x) unlist(strsplit(x,split=".p"))[1])

#Identify results files
results_files <- grep(sleuth_files,pattern="results.tsv",value = TRUE)

#Read in each result file, add to list then bind all into one tibble

all_res_list <- list()
for (i in 1:length(results_files)){
  all_res_list[[i]] <- read_tsv(paste0("08_inputs_transcriptomics/sleuth/",results_files[i]),col_types = c("cddddddcccic"))
}

all_res <- bind_rows(all_res_list)

#Need to translate all gene IDs to human ensembl gene IDs
species_info <- read_csv("08_inputs_transcriptomics/species_info.csv")

ensembl = useMart("ensembl")

translation_list <- list()

#For each species, except for canary and human, get ensembl ids and human ortholog ensembl ids, add each to list
#Note for some reason it says some species dataset names are not valid even though they are. Running again will usualy fix the problem.
for (i in 1:nrow(species_info)){
  species <- species_info[i,"species"]%>%pull
  if (!(species %in% c("Serinus_canaria","Homo_sapiens","Gallus_gallus"))){
    ensembl_sp <- useDataset(species_info[i,"biomart_dataset"]%>%pull,mart=ensembl)
    translation_list[[i]] <- getBM(attributes=c("ensembl_gene_id","hsapiens_homolog_ensembl_gene","ggallus_homolog_ensembl_gene"),
                          #filters=c("with_hsapiens_homolog"),
                          values=TRUE,
                          mart=ensembl_sp) %>% as.tibble
  }
  if (species == "Gallus_gallus"){
    ensembl_sp <- useDataset(species_info[i,"biomart_dataset"]%>%pull,mart=ensembl)
    translation_list[[i]] <- getBM(attributes=c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"),
                                   #filters=c("with_hsapiens_homolog"),
                                   values=TRUE,
                                   mart=ensembl_sp) %>% as.tibble
  }
}

#For canary, need to get ggallus orthologs from hog inputs, then use the previously created gallus gallus (from script 02) translation table to go from chicken NCBI gene IDs to human ensembl IDs

load("02_output_annotated_data/biomart_translation_tables.Rdat")
#Read in HOG IDs
hog_ids <- read_tsv("02_input_annotation_data/new_hog_list.txt",col_names =  c("HOG2_HogID","NCBI_ID","entrezgene","sp"))

#Data cleanup, extract chicken gene IDs
hog_ids_gg <- hog_ids %>%
  separate(HOG2_HogID,sep = "_",into=c("drop","hog")) %>%
  mutate(entrezgene = as.character(entrezgene), sp = as.character(sp)) %>%
  filter(sp == "galGal") %>%
  dplyr::select(hog,entrezgene,sp) %>%
  left_join(gg_trans_table,by=c("entrezgene")) %>%
  dplyr::select(hog,entrezgene,sp,hsapiens_homolog_ensembl_gene,ggallus_homolog_ensembl_gene=ensembl_gene_id)

hog_ids_sc <- hog_ids %>%
  separate(HOG2_HogID,sep = "_",into=c("drop","hog")) %>%
  mutate(entrezgene = as.character(entrezgene), sp = as.character(sp)) %>%
  filter(sp == "serCan") %>%
  dplyr::select(hog,entrezgene_sc = entrezgene,sp)

#Combine and select columns to add to overall translation list
hog_ids_sc_combo <- hog_ids_sc %>%
  left_join(hog_ids_gg,by="hog") %>%
  filter(!is.na(hsapiens_homolog_ensembl_gene)) %>%
  dplyr::select(ensembl_gene_id = entrezgene_sc, hsapiens_homolog_ensembl_gene,ggallus_homolog_ensembl_gene)

#Idenfity correct index in translation list, and add results
index <- grep("Serinus_canaria",species_info$species,value=F)

translation_list[[index]] <- hog_ids_sc_combo

#Fix chicken spot to include gallus ensembl ids for ggallus_homolog_ensembl_gene
index <- grep("Gallus_gallus",species_info$species,value=F)
translation_list[[index]] <- translation_list[[index]] %>%
  mutate(ggallus_homolog_ensembl_gene=ensembl_gene_id)

#Create table to use for all species
trans_table <- bind_rows(translation_list)

write_csv(trans_table,"08_output_transcriptomics/all_species_translation_table.csv")


#Read in infectious agent info
infect_info <- read_csv("08_inputs_transcriptomics/infect_info_table.csv")


#########################################################################################################
################### Bird Only Transcriptome dataset cleanup #############################################
#########################################################################################################
#First, going to analyze bird-only data with large bird selection dataset, to take advantage of all of the bird genes we can't match to humans.

#First, need to get rid of extra info on end of target_ids, and create new column of sig given beta > 1 or -1
all_res_birds <- all_res %>%
  separate(target_id,into=c("ensembl_id","drop"),fill = "right") %>%
  dplyr::select(-drop) %>%
  unite("bioproject_infect",c("bioproject","infect"),sep = "-",remove=FALSE) %>%
  dplyr::select(-bioproject) %>%
  dplyr::rename(bioproject = bioproject_infect) %>%
  left_join(trans_table,by=c("ensembl_id" = "ensembl_gene_id")) %>%
  mutate(ensembl_id_gg = if_else(species=="GAL",ensembl_id,ggallus_homolog_ensembl_gene)) %>%
  mutate(beta_TF = if_else(abs(b)>=1,1,0), beta_sig = sig*beta_TF) %>%
  separate(bioproject, into=c("bioproj_number"), extra="drop",remove = F,sep = "-") %>%
  left_join(infect_info,by=c("infect","species")) %>%
  filter(clade == "birds")

#Now, want to make sure there are not duplicate gene IDs in each bioproject. For those that are, take the sum of beta_sig and sig, and change any >1 to 1 and <-1 to -1.
all_res_birds_nodups <- all_res_birds %>%
  group_by(bioproject,ensembl_id_gg) %>%
  summarize(beta_sig = sum(beta_sig), sig = sum(sig)) %>%
  mutate(singles_beta_sig = case_when(beta_sig >= 1 ~ 1,
                              beta_sig == 0 ~ 0,
                              beta_sig <= -1 ~ -1),
         singles_sig = case_when(sig >= 1 ~ 1,
                         sig == 0 ~ 0,
                         sig <= -1 ~ -1)) %>%
  dplyr::select(-beta_sig,-sig) %>%
  ungroup()

#Combine back with all_res_birds, remove duplicates and drop old columns
all_res_birds_singles <- all_res_birds %>%
  left_join(all_res_birds_nodups,by=c("bioproject","ensembl_id_gg")) %>%
  dplyr::select(-beta_sig,-sig) %>%
  rename(beta_sig = singles_beta_sig,sig=singles_sig) %>%
  distinct(bioproject,ensembl_id_gg,.keep_all=T)

#Now each specific dataset is cleaned up, but we need to clean up bioprojects that had more than 1 condition (timepoints, etc.)
bioprojects_n_conditions <- all_res_birds_singles %>%
  group_by(bioproj_number,bioproject) %>%
  distinct(bioproj) %>%
  group_by(bioproj_number) %>%
  summarize(n_cond = n())

#First, create a clean dataset of those bioprojects that only have one condition:
single_bioprojects <- bioprojects_n_conditions %>%
  filter(n_cond == 1) %>%
  pull(bioproj_number)

all_res_birds_singles_clean <- all_res_birds_singles %>%
  filter(bioproj_number %in% single_bioprojects)



#Now, for those that have only two conditions (e.g. different organs). Want to include genes that are significant in either condition 
double_bioprojects <- bioprojects_n_conditions %>%
  filter(n_cond == 2) %>%
  pull(bioproj_number)

#First check, are there any genes that are upregulated in one condition and down in the other?
all_res_birds_singles %>%
  filter(bioproj_number %in% double_bioprojects) %>%
  group_by(bioproj_number,ensembl_id_gg) %>%
  summarize(beta_sig = sum(beta_sig), sig = sum(sig),beta_sig_abs = sum(abs(beta_sig)),sig_abs = sum(abs(sig))) %>%
  filter(abs(beta_sig) != beta_sig_abs | abs(sig) != sig_abs)

#No, so we don't have to worry about that
sig_values_birds_doubles_fixed <- all_res_birds_singles %>%
  filter(bioproj_number %in% double_bioprojects) %>%
  group_by(bioproj_number,ensembl_id_gg) %>%
  summarize(beta_sig = sum(beta_sig), sig = sum(sig))  %>%
  mutate(singles_beta_sig = case_when(beta_sig >= 1 ~ 1,
                                      beta_sig == 0 ~ 0,
                                      beta_sig <= -1 ~ -1),
         singles_sig = case_when(sig >= 1 ~ 1,
                                 sig == 0 ~ 0,
                                 sig <= -1 ~ -1)) %>%
  dplyr::select(-beta_sig,-sig) %>%
  ungroup()

#Combine back with all_res_birds, remove duplicates and drop old columns
all_res_birds_doubles_clean <-all_res_birds_singles %>%
  filter(bioproj_number %in% double_bioprojects) %>%
  left_join(sig_values_birds_doubles_fixed,by=c("bioproj_number","ensembl_id_gg")) %>%
  dplyr::select(-beta_sig,-sig) %>%
  rename(beta_sig = singles_beta_sig,sig=singles_sig) %>%
  distinct(bioproj_number,ensembl_id_gg,.keep_all=T)



#Next we will fix all of the bioprojects that have 4 different time points/conditions. We will do the same procedure as above, but consider a gene significantly up or down regulated if it is in at least 2 conditions.
quad_bioprojects <- bioprojects_n_conditions %>%
  filter(n_cond == 4) %>%
  pull(bioproj_number)

#First check, are there any genes that are upregulated in one condition and down in the other?
all_res_birds_singles %>%
  filter(bioproj_number %in% quad_bioprojects) %>%
  group_by(bioproj_number,ensembl_id_gg) %>%
  summarize(beta_sig = sum(beta_sig), sig = sum(sig),beta_sig_abs = sum(abs(beta_sig)),sig_abs = sum(abs(sig))) %>%
  filter(abs(beta_sig) != beta_sig_abs | abs(sig) != sig_abs)
#No, so we don't have to worry about that

sig_values_birds_quad_fixed <- all_res_birds_singles %>%
  filter(bioproj_number %in% quad_bioprojects) %>%
  group_by(bioproj_number,ensembl_id_gg) %>%
  summarize(beta_sig = sum(beta_sig), sig = sum(sig))  %>%
  mutate(singles_beta_sig = case_when(beta_sig >= 2 ~ 1,
                                      beta_sig < 2 & beta_sig >-2 ~ 0,
                                      beta_sig <= -2 ~ -1),
         singles_sig = case_when(sig >= 2 ~ 1,
                                 sig < 2 & sig >-2 ~ 0,
                                 sig <= -2 ~ -1)) %>%
  dplyr::select(-beta_sig,-sig) %>%
  ungroup()

#Combine back with all_res_birds, remove duplicates and drop old columns
all_res_birds_quad_clean <-all_res_birds_singles %>%
  filter(bioproj_number %in% quad_bioprojects) %>%
  left_join(sig_values_birds_quad_fixed,by=c("bioproj_number","ensembl_id_gg")) %>%
  dplyr::select(-beta_sig,-sig) %>%
  rename(beta_sig = singles_beta_sig,sig=singles_sig) %>%
  distinct(bioproj_number,ensembl_id_gg,.keep_all=T)



#Next we will fix all of the bioprojects that have 6 different time points/conditions. We will do the same procedure as above, but consider a gene significantly up or down regulated if it is in at least 3 conditions.
six_bioprojects <- bioprojects_n_conditions %>%
  filter(n_cond == 6) %>%
  pull(bioproj_number)

#First check, are there any genes that are upregulated in one condition and down in the other?
all_res_birds_singles %>%
  filter(bioproj_number %in% six_bioprojects) %>%
  group_by(bioproj_number,ensembl_id_gg) %>%
  summarize(beta_sig = sum(beta_sig), sig = sum(sig),beta_sig_abs = sum(abs(beta_sig)),sig_abs = sum(abs(sig))) %>%
  filter(abs(beta_sig) != beta_sig_abs | abs(sig) != sig_abs)
#No, so we don't have to worry about that

sig_values_birds_six_fixed <- all_res_birds_singles %>%
  filter(bioproj_number %in% six_bioprojects) %>%
  group_by(bioproj_number,ensembl_id_gg) %>%
  summarize(beta_sig = sum(beta_sig), sig = sum(sig))  %>%
  mutate(singles_beta_sig = case_when(beta_sig >= 3 ~ 1,
                                      beta_sig < 3 & beta_sig >-3 ~ 0,
                                      beta_sig <= -3 ~ -1),
         singles_sig = case_when(sig >= 3 ~ 1,
                                 sig < 3 & sig >-3 ~ 0,
                                 sig <= -3 ~ -1)) %>%
  dplyr::select(-beta_sig,-sig) %>%
  ungroup()

#Combine back with all_res_birds, remove duplicates and drop old columns
all_res_birds_six_clean <-all_res_birds_singles %>%
  filter(bioproj_number %in% six_bioprojects) %>%
  left_join(sig_values_birds_six_fixed,by=c("bioproj_number","ensembl_id_gg")) %>%
  dplyr::select(-beta_sig,-sig) %>%
  rename(beta_sig = singles_beta_sig,sig=singles_sig) %>%
  distinct(bioproj_number,ensembl_id_gg,.keep_all=T)



#Next we will fix all of the bioprojects that have 8 different time points/conditions. We will do the same procedure as above, but consider a gene significantly up or down regulated if it is in at least 2 conditions. Here I am being more lenient, because many of the treatments had very few significantly differently expressed genes.
eight_bioprojects <- bioprojects_n_conditions %>%
  filter(n_cond == 8) %>%
  pull(bioproj_number)

#First check, are there any genes that are upregulated in one condition and down in the other?
all_res_birds_singles %>%
  filter(bioproj_number %in% eight_bioprojects) %>%
  group_by(bioproj_number,ensembl_id_gg) %>%
  summarize(beta_sig = sum(beta_sig), sig = sum(sig),beta_sig_abs = sum(abs(beta_sig)),sig_abs = sum(abs(sig))) %>%
  filter(abs(beta_sig) != beta_sig_abs | abs(sig) != sig_abs)
#No, so we don't have to worry about that

sig_values_birds_eight_fixed <- all_res_birds_singles %>%
  filter(bioproj_number %in% eight_bioprojects) %>%
  group_by(bioproj_number,ensembl_id_gg) %>%
  summarize(beta_sig = sum(beta_sig), sig = sum(sig)) %>%
  mutate(singles_beta_sig = case_when(beta_sig >= 2 ~ 1,
                                      beta_sig < 2 & beta_sig >=-1 ~ 0,
                                      (beta_sig <= -2) ~ -1),
         singles_sig = case_when(sig >= 2 ~ 1,
                                 sig < 2 & sig >-2 ~ 0,
                                 sig <= -2 ~ -1)) %>%
  dplyr::select(-beta_sig,-sig) %>%
  ungroup()

#Combine back with all_res_birds, remove duplicates and drop old columns
all_res_birds_eight_clean <-all_res_birds_singles %>%
  filter(bioproj_number %in% eight_bioprojects) %>%
  left_join(sig_values_birds_eight_fixed,by=c("bioproj_number","ensembl_id_gg")) %>%
  dplyr::select(-beta_sig,-sig) %>%
  rename(beta_sig = singles_beta_sig,sig=singles_sig) %>%
  distinct(bioproj_number,ensembl_id_gg,.keep_all=T)






#Combine all results into a single, clean tibble
all_res_birds_clean <- bind_rows(all_res_birds_singles_clean,
                                  all_res_birds_doubles_clean,
                                  all_res_birds_quad_clean,
                                  all_res_birds_six_clean,
                                  all_res_birds_eight_clean) %>%
  dplyr::select(-bioproject)

save(all_res_birds_clean,file = "08_output_transcriptomics/birds_clean_transcriptomic_results.Rdat")


#What are the bioprojects by infections?
bioproj_infect_info <- all_res_birds_clean %>%
  distinct(bioproj_number,.keep_all=T) %>%
  dplyr::select(bioproj_number,infect_type_specific,infect_type_virus_type,infect_type_general)

#Load in the significance results
load("02_output_annotated_data/all_res_zf_hs.Rdat")

#First combine the transcriptomic results with the gene tree significance results, cut down some unnecessary columns in the process, note that our dataset shrinks a bit because not all hogs could be assigned gallus gallus ensembl gene ids. Filled in those with the zf ensembl id. Numbers match those that have NCBI gene IDs
sig_res_simple <- all_res_gene_zf_hs %>%
  dplyr::select(hog,ensembl_id_gg = ensembl_gene_id,ensembl_id_zf=ensembl_gene_id_zf,FDRPval_m1m2,FDRPval_m2m2a,FDRPval_m7m8,FDRPval_m8m8a,FDRPval_busted,prop_sel.n,external_gene_name=external_gene_name.x,entrezgene) %>%
  left_join(trans_table,by=c("ensembl_id_zf"="ensembl_gene_id")) %>%
  mutate(ensembl_id_gg=if_else(is.na(ensembl_id_gg),ggallus_homolog_ensembl_gene,ensembl_id_gg)) %>%
  distinct(ensembl_id_gg,.keep_all=TRUE) %>%
  filter(!is.na(ensembl_id_gg)) %>%
  mutate(sig_all_tests = if_else(FDRPval_m1m2 < 0.05 & FDRPval_m2m2a < 0.05 & FDRPval_m7m8 < 0.05 & FDRPval_m8m8a < 0.05 & FDRPval_busted < 0.05,1,0)) %>%
  dplyr::select(ensembl_id_gg,prop_sel.n,external_gene_name,sig_all_tests)

#Combine those significance results with transcriptome results

infect_agents <- all_res_birds_clean %>%
  distinct(infect_type_virus_type) %>%
  pull

#Colors to use:
general_colors_plasmodium <- c("white",rgb(27,158,119,225,maxColorValue = 225))
names(general_colors_plasmodium) <- c("not significant","significant")

general_colors_virus <- c("white",rgb(117,112,179,225,maxColorValue = 225))
names(general_colors_virus) <- c("not significant","significant")

general_colors_bacterium<- c("white",rgb(217,95,2,225,maxColorValue = 225))
names(general_colors_bacterium) <- c("not significant","significant")

color_list <- list(general_colors_plasmodium,general_colors_virus,general_colors_bacterium)
names(color_list) <- c("plasmodium","virus","bacterium")

color_vec <- as.tibble(color_list)[2,] %>% gather(agent,color) %>% pull(color)
names(color_vec) <- as.tibble(color_list)[2,] %>% gather(agent,color) %>% pull(agent)



infect_res_table <- matrix(nrow=length(infect_agents)*3,ncol=8)
colnames(infect_res_table) <- c("infect_agent", "trans_response", "n", "estimate","std.error","z.value","p.value", "prop_sig")

infect_plots <- list()

start <- 1

for (i in 1:length(infect_agents)){
  infect_res_table[(start:(start+2)),1] <- infect_agents[i]
  infect_res <- all_res_birds_clean %>%
    filter(infect_type_virus_type==infect_agents[i]) %>%
    group_by(ensembl_id_gg) %>%
    summarize(beta_sig_overall = sum(beta_sig), sig_overall = sum(sig)) %>%
    mutate(beta_sig = case_when(beta_sig_overall >= 1 ~ 1,
                                beta_sig_overall == 0 ~ 0,
                                beta_sig_overall <= -1 ~ -1),
           sig = case_when(sig_overall >= 1 ~ 1,
                           sig_overall == 0 ~ 0,
                           sig_overall <= -1 ~ -1)) %>%
    filter(ensembl_id_gg != "", !is.na(beta_sig)) %>%
    left_join(sig_res_simple,by=c("ensembl_id_gg")) %>%
    filter(!is.na(sig_all_tests)) %>%
    mutate(expr_sig = if_else(beta_sig !=0,1,0), up_reg = if_else(beta_sig==1,1,0), down_reg = if_else(beta_sig==-1,1,0), not_reg = if_else(beta_sig == 0,1,0))
  
  infect_res_table[start,2] <- "down"   
  try(infect_res_downreg_test <- infect_res %>%
    filter(beta_sig != 1) %>%
    glm(sig_all_tests ~ down_reg, family="binomial",data=.))
  try(infect_res_table[start,4:7] <- summary(infect_res_downreg_test)$coefficient[2,])
  
  infect_res_table[start+1,2] <- "none" 
  infect_res_table[start+1,4:7] <- NA

  infect_res_table[start+2,2] <- "up"   
  try(infect_res_upreg_test <- infect_res %>%
    filter(beta_sig != -1) %>%
    glm(sig_all_tests ~ up_reg, family="binomial",data=.))
  try(infect_res_table[start+2,4:7] <- summary(infect_res_upreg_test)$coefficient[2,])
  
  #infect_res %>%
  #  filter(up_reg ==1, sig_all_tests ==1) %>%
  #  print(n=30)
  
  general_type <- bioproj_infect_info %>%
    filter(infect_type_virus_type == infect_agents[i]) %>%
    distinct(infect_type_general) %>%
    pull
  
  #Calculate proportion significant
  infect_res_table[start:(start+2),8] <- infect_res %>%
    with(.,table(sig_all_tests,beta_sig)) %>%
    as.tibble %>%
    mutate(beta_sig = case_when(beta_sig == -1 ~ "down",
                                beta_sig == 0 ~ "none",
                                beta_sig == 1 ~ "up")) %>%
    group_by(beta_sig) %>%
    mutate(n_cat = sum(n)) %>%
    spread(sig_all_tests,n) %>%
    mutate(proportion = `1`/n_cat) %>%
    pull(proportion)
  
  #Pull the number of genes in each transcriptional category
  infect_res_table[start:(start+2),3] <- infect_res %>%
    with(.,table(beta_sig)) %>%
    as.tibble %>%
    mutate(beta_sig = case_when(beta_sig == -1 ~ "down",
                                beta_sig == 0 ~ "none",
                                beta_sig == 1 ~ "up")) %>%
    pull(n)
  
  
  
  infect_plots[[i]] <-  infect_res_table[start:(start+2),] %>%
    as.tibble %>%
    mutate_at(.vars = 3:8, .fun = as.numeric) %>%
    mutate(sig_char = case_when(p.value <= 0.0001 ~ "***",
                                p.value > 0.0001 & p.value <= 0.001 ~ "**",
                                p.value > 0.001 & p.value <= 0.05 ~ "*",
                                p.value > 0.05 ~ "")) %>%
    ggplot(aes(trans_response,prop_sig)) +
    geom_bar(stat="identity",fill=color_vec[general_type]) +
    ggtitle(infect_agents[i]) +
    ylab("proportion significant") +
    xlab("transcriptional response") +
    ylim(0,1) +
    geom_text(aes(trans_response,label=n),nudge_y=0.08,size=4) +
    geom_text(aes(trans_response,label=sig_char),nudge_y=0.01,size=4)

  start <- start + 3 
  }

names(infect_plots) <- infect_agents

birds_legend <- all_res_birds_clean %>%
  ggplot(aes(infect_type_general,fill=factor(infect_type_general,levels=c("plasmodium","virus","bacterium")))) +
  geom_bar(position="fill") +
  scale_fill_manual(values = color_vec,name="pathogen type") +
  ylab("proportion") +
  ylim(1,2) +
  guides(fill=guide_legend(override.aes = list(fill=color_vec))) +
  theme(line = element_blank(),
        axis.text = element_blank(),
        title = element_blank())

plot_grid(infect_plots[["birnaviridae"]],infect_plots[["influenza"]],infect_plots[["paramyxovirus"]],infect_plots[["west_nile_virus"]],infect_plots[["ecoli"]],infect_plots[["mycoplasma"]],infect_plots[["plasmodium"]],birds_legend,ncol=4)
ggsave("08_output_transcriptomics/bird_transcriptome_sig_figure.pdf",height=7,width=11)








################ Need to translate to human gene IDs 
#First, need to get rid of extra info on end of target_ids, and create new column of sig given beta > 1 or -1

all_res_anno <- all_res %>%
  separate(target_id,into=c("ensembl_id","drop"),fill = "right") %>%
  dplyr::select(-drop) %>%
  left_join(trans_table,by=c("ensembl_id" = "ensembl_gene_id")) %>%
  mutate(ensembl_id_hs = if_else(species=="HOM",ensembl_id,hsapiens_homolog_ensembl_gene)) %>%
  mutate(beta_TF = if_else(abs(b)>=1,1,0), beta_sig = sig*beta_TF)

save(all_res_anno,file = "08_output_transcriptomics/all_transcriptomic_results_human_ensembl_ids.Rdat")

#How many ensembl_ids don't have a human match?
n_matches <- all_res_anno %>%
  group_by(bioproject,species) %>%
  summarize(n_transcripts = n(),n_no_match = sum(is.na(ensembl_id_hs)))

#Further split up bioproject for easier analysis
all_res_anno <- all_res_anno %>%
  separate(bioproject, into=c("bioproj_number"), extra="drop",remove = F,sep = "-") %>%
  left_join(infect_info,by=c("infect","species"))

#Create a tibble with useful info for each bioproject
bioproject_info <- all_res_anno %>%
  distinct(bioproject,.keep_all=T) %>%
  dplyr::select(bioproject,bioproj_number,infect,infect_type_specific,infect_type_virus_type,infect_type_general,clade)

all_res_anno %>%
  filter(!is.na(sig)) %>%
  group_by(bioproject) %>%
  summarize(n_up = sum(sig==1),n_up_beta = sum(beta_sig==1), n_down = sum(sig==-1), n_down_beta = sum(beta_sig==-1), n_notsig = sum(sig==0), n_notsig_beta = sum(beta_sig ==0),n_ensID_hs = sum(!is.na(ensembl_id_hs))) %>%
  left_join(bioproject_info, by="bioproject") %>%
  write_csv("08_output_transcriptomics/bioproject_basic_info.csv")

#Now, to clean up each bioproject. 
#For each dataset (bioproject), make sure each geneID is only represented once. Will sum across gene IDs, and turn any >1 into 1, and any <-1 into -1.
all_res_anno %>%
  group_by(bioproject, ensembl_id_hs) %>%
  


#First, bioprojects with single result files (no timepoints, differences in infectious agents, etc.)



#Read in significance results
imm <- read_csv("05_output_bird_mammal_comparison_results/bird_mammal_combined_dataset.csv")

qvals <- c(0.1,0.01,0.001,0.0001)
qval <- 0.01












##########Test first for infection_type_general, then infect_type_virus_type. Loop through and calc similar up-regulated, down-regulated, and just significant. Test for relationships (logistic regression) for significant genes in birds given sig genes in mammals and transcriptome status, and vica versa for birds and mammals.
general_types <- all_res_anno %>%
  distinct(infect_type_general) %>%
  pull

trans_tests <- c("sig_notsig","up_reg","down_reg")

mammals_general_types_res <- matrix(nrow=length(general_types)*length(trans_tests)*8,ncol=8)
birds_general_types_res <- matrix(nrow=length(general_types)*length(trans_tests)*8,ncol=8)


start <- 1

for (i in 1:length(general_types)){
all_res_anno_test <- all_res_anno %>%
  filter(infect_type_general == general_types[i])

  for (j in 1:length(trans_tests)){
    mammals_general_types_res[start:(start+7),1] <- general_types[i]
    birds_general_types_res[start:(start+7),1] <- general_types[i]
    mammals_general_types_res[start:(start+7),2] <- trans_tests[j]
    birds_general_types_res[start:(start+7),2] <- trans_tests[j]
  #Get sig counts for each gene for both birds and mammals
    ens_sig_birds_mammals_test <- all_res_anno_test %>%
      filter(!is.na(ensembl_id_hs)) %>%
      mutate(sig_notsig = if_else(sig != 0, 1,0),up_reg = if_else(sig==1,1,0),down_reg = if_else(sig ==-1,1,0)) %>%
      group_by(ensembl_id_hs,clade) %>%
      summarize(n_sig_notsig = sum(sig_notsig), n_up_reg = sum(up_reg), n_down_reg = sum(down_reg)) %>%
      dplyr::select(ensembl_id_hs,clade,paste0("n_",trans_tests[j])) %>%
      spread(clade,paste0("n_",trans_tests[j]))


#Filter out missing genes in birds, make a table with counts in each category to get a sense of the data
    imm_test <- imm %>%
      left_join(ens_sig_birds_mammals_test, by=c("ensembl_gene_id_hs"="ensembl_id_hs")) %>%
      filter(!is.na(birds), !is.na(mammals)) %>%
      mutate(status = case_when(
        birds>0 & mammals==0 ~ "birds_only",
        birds==0 & mammals>0 ~ "mammals_only",
        birds>0 & mammals>0 ~ "birds_and_mammals",
        birds==0 & mammals==0 ~ "neither"),
        status_mammals = if_else(mammals>0,TRUE,FALSE),
        status_birds = if_else(birds>0,TRUE,FALSE))
    imm_test %>%
      with(.,table(status))

    gl_res_mammals <- glm(sig_mammals ~ sig_birds*status, family="binomial",data=imm_test)
    #summary(gl_res_mammals)
    mammals_general_types_res[start:(start+6),3] <- "sig_mammals"
    mammals_general_types_res[start:(start+6),4] <- rownames(summary(gl_res_mammals)$coefficient[2:8,])
    mammals_general_types_res[start:(start+6),5:8] <- summary(gl_res_mammals)$coefficient[2:8,]
    
    #Test if sig transcriptome results correspond to sig selection results
    gl_res_mammals_only <- glm(sig_mammals ~ status_mammals, family="binomial",data=imm_test)
    mammals_general_types_res[(start+7),4] <- rownames(summary(gl_res_mammals_only)$coefficient)[2]
    mammals_general_types_res[(start+7),5:8] <- summary(gl_res_mammals_only)$coefficient[2,]
    

    gl_res_birds <- glm(sig_birds ~ sig_mammals*status, family="binomial",data=imm_test)
    #summary(gl_res_birds)
    birds_general_types_res[start:(start+6),3] <- "sig_birds"
    birds_general_types_res[start:(start+6),4] <- rownames(summary(gl_res_birds)$coefficient[2:8,])
    birds_general_types_res[start:(start+6),5:8] <- summary(gl_res_birds)$coefficient[2:8,]
    #Test if sig transcriptome results correspond to sig selection results
    gl_res_birds_only <- glm(sig_birds ~ status_birds, family="binomial",data=imm_test)
    birds_general_types_res[(start+7),4] <- rownames(summary(gl_res_birds_only)$coefficient)[2]
    birds_general_types_res[(start+7),5:8] <- summary(gl_res_birds_only)$coefficient[2,]
    

    start <- start + 8

    }
  }

colnames(mammals_general_types_res) <- c("agent_class","trans_comparison","dep_var","parameter","estimate","std.error","z.value","p.value")
colnames(birds_general_types_res) <- c("agent_class","trans_comparison","dep_var","parameter","estimate","std.error","z.value","p.value")

#Clean up, select relevant columns
mammals_general_types_res_clean <- mammals_general_types_res %>%
  as.tibble %>%
  mutate(p.value=round(as.numeric(p.value),digits = 4),
         estimate=round(as.double(estimate),digits=2),
         std.error=round(as.numeric(std.error),digits=2),
         z.value=round(as.numeric(z.value),digits=2)) %>%
  dplyr::select(agent_class,trans_comparison,parameter,estimate,std.error,z.value,p.value)
birds_general_types_res_clean <- birds_general_types_res %>%
  as.tibble %>%
  mutate(p.value=round(as.numeric(p.value),digits = 4),
         estimate=round(as.double(estimate),digits=2),
         std.error=round(as.numeric(std.error),digits=2),
         z.value=round(as.numeric(z.value),digits=2)) %>%
  dplyr::select(agent_class,trans_comparison,parameter,estimate,std.error,z.value,p.value)

write_csv(mammals_general_types_res_clean,"08_output_transcriptomics/general_class_mammal_logistic_regressions.csv")
write_csv(birds_general_types_res_clean,"08_output_transcriptomics/general_class_bird_logistic_regressions.csv")









##########Next, test for infect_type_virus_type. Loop through and calc similar up-regulated, down-regulated, and just significant. Test for relationships (logistic regression) for significant genes in birds given sig genes in mammals and transcriptome status, and vica versa for birds and mammals.
virus_types <- all_res_anno %>%
  distinct(infect_type_virus_type) %>%
  pull

trans_tests <- c("sig_notsig","up_reg","down_reg")

mammals_virus_types_res <- matrix(nrow=length(virus_types)*length(trans_tests)*8,ncol=8)
birds_virus_types_res <- matrix(nrow=length(virus_types)*length(trans_tests)*8,ncol=8)


start <- 1

for (i in 1:length(virus_types)){
  all_res_anno_test <- all_res_anno %>%
    filter(infect_type_virus_type == virus_types[i])
  
  for (j in 1:length(trans_tests)){
    mammals_virus_types_res[start:(start+7),1] <- virus_types[i]
    birds_virus_types_res[start:(start+7),1] <- virus_types[i]
    mammals_virus_types_res[start:(start+7),2] <- trans_tests[j]
    birds_virus_types_res[start:(start+7),2] <- trans_tests[j]
    #Get sig counts for each gene for both birds and mammals
    ens_sig_birds_mammals_test <- all_res_anno_test %>%
      filter(!is.na(ensembl_id_hs)) %>%
      mutate(sig_notsig = if_else(sig != 0, 1,0),up_reg = if_else(sig==1,1,0),down_reg = if_else(sig ==-1,1,0)) %>%
      group_by(ensembl_id_hs,clade) %>%
      summarize(n_sig_notsig = sum(sig_notsig), n_up_reg = sum(up_reg), n_down_reg = sum(down_reg)) %>%
      dplyr::select(ensembl_id_hs,clade,paste0("n_",trans_tests[j])) %>%
      spread(clade,paste0("n_",trans_tests[j]))
    
    
    #Filter out missing genes in birds, make a table with counts in each category to get a sense of the data
    imm_test <- imm %>%
      left_join(ens_sig_birds_mammals_test, by=c("ensembl_gene_id_hs"="ensembl_id_hs")) %>%
      filter(!is.na(birds), !is.na(mammals)) %>%
      mutate(status = case_when(
        birds>0 & mammals==0 ~ "birds_only",
        birds==0 & mammals>0 ~ "mammals_only",
        birds>0 & mammals>0 ~ "birds_and_mammals",
        birds==0 & mammals==0 ~ "neither"),
        status_mammals = if_else(mammals>0,TRUE,FALSE),
        status_birds = if_else(birds>0,TRUE,FALSE))
    imm_test %>%
      with(.,table(status))
    
    gl_res_mammals <- glm(sig_mammals ~ sig_birds*status, family="binomial",data=imm_test)
    #summary(gl_res_mammals)
    mammals_virus_types_res[start:(start+6),3] <- "sig_mammals"
    try(mammals_virus_types_res[start:(start+6),4] <- rownames(summary(gl_res_mammals)$coefficient[2:8,]))
    try(mammals_virus_types_res[start:(start+6),5:8] <- summary(gl_res_mammals)$coefficient[2:8,])
    
    #Test if sig transcriptome results correspond to sig selection results
    gl_res_mammals_only <- glm(sig_mammals ~ status_mammals, family="binomial",data=imm_test)
    try(mammals_virus_types_res[(start+7),4] <- rownames(summary(gl_res_mammals_only)$coefficient)[2])
    try(mammals_virus_types_res[(start+7),5:8] <- summary(gl_res_mammals_only)$coefficient[2,])
    
    
    gl_res_birds <- glm(sig_birds ~ sig_mammals*status, family="binomial",data=imm_test)
    #summary(gl_res_birds)
    birds_virus_types_res[start:(start+6),3] <- "sig_birds"
    try(birds_virus_types_res[start:(start+6),4] <- rownames(summary(gl_res_birds)$coefficient[2:8,]))
    try(birds_virus_types_res[start:(start+6),5:8] <- summary(gl_res_birds)$coefficient[2:8,])
    #Test if sig transcriptome results correspond to sig selection results
    gl_res_birds_only <- glm(sig_birds ~ status_birds, family="binomial",data=imm_test)
    try(birds_virus_types_res[(start+7),4] <- rownames(summary(gl_res_birds_only)$coefficient)[2])
    try(birds_virus_types_res[(start+7),5:8] <- summary(gl_res_birds_only)$coefficient[2,])
    
    
    start <- start + 8
    
  }
}

colnames(mammals_virus_types_res) <- c("agent_class","trans_comparison","dep_var","parameter","estimate","std.error","z.value","p.value")
colnames(birds_virus_types_res) <- c("agent_class","trans_comparison","dep_var","parameter","estimate","std.error","z.value","p.value")

#Clean up, select relevant columns
mammals_virus_types_res_clean <- mammals_virus_types_res %>%
  as.tibble %>%
  mutate(p.value=round(as.numeric(p.value),digits = 4),
         estimate=round(as.double(estimate),digits=2),
         std.error=round(as.numeric(std.error),digits=2),
         z.value=round(as.numeric(z.value),digits=2)) %>%
  dplyr::select(agent_class,trans_comparison,parameter,estimate,std.error,z.value,p.value)
birds_virus_types_res_clean <- birds_virus_types_res %>%
  as.tibble %>%
  mutate(p.value=round(as.numeric(p.value),digits = 4),
         estimate=round(as.double(estimate),digits=2),
         std.error=round(as.numeric(std.error),digits=2),
         z.value=round(as.numeric(z.value),digits=2)) %>%
  dplyr::select(agent_class,trans_comparison,parameter,estimate,std.error,z.value,p.value)

write_csv(mammals_virus_types_res_clean,"08_output_transcriptomics/virus_class_mammal_logistic_regressions.csv")
write_csv(birds_virus_types_res_clean,"08_output_transcriptomics/virus_class_bird_logistic_regressions.csv")






























#First look at plasmodium - this is simple, no time points, just a single bioproject for each
all_res_anno_plas <- all_res_anno %>%
  filter(infect_type_general == "plasmodium")

#Get sig counts for each gene for both birds and mammals
ens_sig_birds_mammals_plas <- all_res_anno_plas %>%
  filter(!is.na(ensembl_id_hs)) %>%
  mutate(sig_notsig = if_else(sig != 0, 1,0),up_reg = if_else(sig==1,1,0),down_reg = if_else(sig ==-1,1,0)) %>%
  group_by(ensembl_id_hs,clade) %>%
  summarize(n_sig = sum(up_reg)) %>%
  spread(clade,n_sig)

#Filter out missing genes in birds, make a table with counts in each category to get a sense of the data
ens_sig_birds_mammals_plas %>%
  filter(!is.na(birds), !is.na(mammals)) %>%
  mutate(status = case_when(
    birds>0 & mammals==0 ~ "birds_only",
    birds==0 & mammals>0 ~ "mammals_only",
    birds>0 & mammals>0 ~ "birds_and_mammals",
    birds==0 & mammals==0 ~ "neither")) %>%
  with(.,table(status))

imm_plas <- imm %>%
  left_join(ens_sig_birds_mammals_plas, by=c("ensembl_gene_id_hs"="ensembl_id_hs")) %>%
  filter(!is.na(birds), !is.na(mammals)) %>%
  mutate(status = case_when(
    birds>0 & mammals==0 ~ "birds_only",
    birds==0 & mammals>0 ~ "mammals_only",
    birds>0 & mammals>0 ~ "birds_and_mammals",
    birds==0 & mammals==0 ~ "neither")) 

imm_plas %>%
  ggplot(aes(status,fill=sig_birds_mammals)) +
  geom_bar(position="fill")
imm_plas %>%
  mutate(birds_reg = if_else(birds > 0,TRUE,FALSE)) %>%
  ggplot(aes(sig_birds,fill=birds_reg)) +
  geom_bar(position="fill")
ggsave("08_output_transcriptomics/birds_plasmodium_sig_regulated.png",height=6,width=4)
imm_plas %>%
  mutate(mammals_reg = if_else(mammals > 0,TRUE,FALSE)) %>%
  ggplot(aes(sig_mammals,fill=mammals_reg)) +
  geom_bar(position="fill")
ggsave("08_output_transcriptomics/mammals_plasmodium_down_regulated.png",height=6,width=4)
ggsave("08_output_transcriptomics/mammals_plasmodium_up_regulated.png",height=6,width=4)
ggsave("08_output_transcriptomics/mammals_plasmodium_sig_regulated.png",height=6,width=4)


#######################Viruses
#Next look at viruses- this is simple, no time points, just a single bioproject for each
all_res_anno_virus <- all_res_anno %>%
  filter(infect_type_general == "virus")

#Get sig counts for each gene for both birds and mammals
ens_sig_birds_mammals_virus <- all_res_anno_virus %>%
  filter(!is.na(ensembl_id_hs)) %>%
  mutate(sig_notsig = if_else(sig != 0, 1,0),up_reg = if_else(sig==1,1,0),down_reg = if_else(sig ==-1,1,0)) %>%
  group_by(ensembl_id_hs,clade) %>%
  summarize(n_sig = sum(down_reg)) %>%
  spread(clade,n_sig)


#Filter out missing genes in birds, make a table with counts in each category to get a sense of the data
imm_virus <- imm %>%
  left_join(ens_sig_birds_mammals_virus, by=c("ensembl_gene_id_hs"="ensembl_id_hs")) %>%
  filter(!is.na(birds), !is.na(mammals)) %>%
  mutate(status = case_when(
    birds>0 & mammals==0 ~ "birds_only",
    birds==0 & mammals>0 ~ "mammals_only",
    birds>0 & mammals>0 ~ "birds_and_mammals",
    birds==0 & mammals==0 ~ "neither")) 
imm_virus %>%
  with(.,table(status))

imm_virus %>%
  ggplot(aes(status,fill=sig_birds_mammals)) +
  geom_bar(position="fill")
imm_virus %>%
  mutate(birds_reg = if_else(birds > 0,TRUE,FALSE)) %>%
  ggplot(aes(sig_birds,fill=birds_reg)) +
  geom_bar(position="fill")
ggsave("08_output_transcriptomics/birds_virus_sig_regulated.png",height=6,width=4)
ggsave("08_output_transcriptomics/birds_virus_up_regulated.png",height=6,width=4)
ggsave("08_output_transcriptomics/birds_virus_down_regulated.png",height=6,width=4)
imm_virus %>%
  mutate(mammals_reg = if_else(mammals > 0,TRUE,FALSE)) %>%
  ggplot(aes(sig_mammals,fill=mammals_reg)) +
  geom_bar(position="fill")
ggsave("08_output_transcriptomics/mammals_virus_sig_regulated.png",height=6,width=4)
ggsave("08_output_transcriptomics/mammals_virus_up_regulated.png",height=6,width=4)
ggsave("08_output_transcriptomics/mammals_virus_down_regulated.png",height=6,width=4)





########Bacteria
#next bacteria- this is simple, no time points, just a single bioproject for each
all_res_anno_bacteria <- all_res_anno %>%
  filter(infect_type_general == "bacterium")

#Get sig counts for each gene for both birds and mammals
ens_sig_birds_mammals_bacteria <- all_res_anno_bacteria %>%
  filter(!is.na(ensembl_id_hs)) %>%
  mutate(sig_notsig = if_else(sig != 0, 1,0),up_reg = if_else(sig==1,1,0),down_reg = if_else(sig ==-1,1,0)) %>%
  group_by(ensembl_id_hs,clade) %>%
  summarize(n_sig = sum(sig_notsig)) %>%
  spread(clade,n_sig)


#Filter out missing genes in birds, make a table with counts in each category to get a sense of the data
imm_bacteria <- imm %>%
  left_join(ens_sig_birds_mammals_bacteria, by=c("ensembl_gene_id_hs"="ensembl_id_hs")) %>%
  filter(!is.na(birds), !is.na(mammals)) %>%
  mutate(status = case_when(
    birds>0 & mammals==0 ~ "birds_only",
    birds==0 & mammals>0 ~ "mammals_only",
    birds>0 & mammals>0 ~ "birds_and_mammals",
    birds==0 & mammals==0 ~ "neither")) 
imm_bacteria %>%
  with(.,table(status))

imm_bacteria %>%
  ggplot(aes(status,fill=sig_birds_mammals)) +
  geom_bar(position="fill")
imm_bacteria %>%
  mutate(birds_reg = if_else(birds > 0,TRUE,FALSE)) %>%
  ggplot(aes(sig_birds,fill=birds_reg)) +
  geom_bar(position="fill")
ggsave("08_output_transcriptomics/birds_bacteria_sig_regulated.png",height=6,width=4)
ggsave("08_output_transcriptomics/birds_bacteria_up_regulated.png",height=6,width=4)
ggsave("08_output_transcriptomics/birds_bacteria_down_regulated.png",height=6,width=4)
imm_bacteria %>%
  mutate(mammals_reg = if_else(mammals > 0,TRUE,FALSE)) %>%
  ggplot(aes(sig_mammals,fill=mammals_reg)) +
  geom_bar(position="fill")
ggsave("08_output_transcriptomics/mammals_bacteria_sig_regulated.png",height=6,width=4)
ggsave("08_output_transcriptomics/mammals_bacteria_up_regulated.png",height=6,width=4)
ggsave("08_output_transcriptomics/mammals_bacteria_down_regulated.png",height=6,width=4)







########################Influenza
#First look at plasmodium - this is simple, no time points, just a single bioproject for each
all_res_anno_influenza <- all_res_anno %>%
  filter(infect_type_virus_type == "influenza")

#Get sig counts for each gene for both birds and mammals
ens_sig_birds_mammals_influenza <- all_res_anno_influenza %>%
  filter(!is.na(ensembl_id_hs)) %>%
  mutate(sig_notsig = if_else(sig != 0, 1,0),up_reg = if_else(sig==1,1,0),down_reg = if_else(sig ==-1,1,0)) %>%
  group_by(ensembl_id_hs,clade) %>%
  summarize(n_sig = sum(sig_notsig)) %>%
  spread(clade,n_sig)


#Filter out missing genes in birds, make a table with counts in each category to get a sense of the data
imm_influenza <- imm %>%
  left_join(ens_sig_birds_mammals_influenza, by=c("ensembl_gene_id_hs"="ensembl_id_hs")) %>%
  filter(!is.na(birds), !is.na(mammals)) %>%
  mutate(status = case_when(
    birds>0 & mammals==0 ~ "birds_only",
    birds==0 & mammals>0 ~ "mammals_only",
    birds>0 & mammals>0 ~ "birds_and_mammals",
    birds==0 & mammals==0 ~ "neither")) 
imm_influenza %>%
  with(.,table(status))

imm_influenza %>%
  ggplot(aes(status,fill=sig_birds_mammals)) +
  geom_bar(position="fill")
imm_influenza %>%
  mutate(birds_reg = if_else(birds > 0,TRUE,FALSE)) %>%
  ggplot(aes(sig_birds,fill=birds_reg)) +
  geom_bar(position="fill")
imm_influenza %>%
  mutate(mammals_reg = if_else(mammals > 0,TRUE,FALSE)) %>%
  ggplot(aes(sig_mammals,fill=mammals_reg)) +
  geom_bar(position="fill")

gl_res_mammals <- glm(sig_mammals ~ sig_birds*status, family="binomial",data=imm_influenza)
summary(gl_res_mammals)

gl_res_birds <- glm(sig_birds ~ sig_mammals*status, family="binomial",data=imm_influenza)
summary(gl_res_birds)




###################Mycoplasma
all_res_anno_mycoplasma <- all_res_anno %>%
  filter(infect_type_virus_type == "mycoplasma")

#Get sig counts for each gene for both birds and mammals
ens_sig_birds_mammals_mycoplasma <- all_res_anno_mycoplasma %>%
  filter(!is.na(ensembl_id_hs)) %>%
  mutate(sig_notsig = if_else(sig != 0, 1,0),up_reg = if_else(sig==1,1,0),down_reg = if_else(sig ==-1,1,0)) %>%
  group_by(ensembl_id_hs,clade) %>%
  summarize(n_sig = sum(sig_notsig)) %>%
  spread(clade,n_sig)

#Filter out missing genes in birds, make a table with counts in each category to get a sense of the data
ens_sig_birds_mammals_mycoplasma %>%
  filter(!is.na(birds), !is.na(mammals)) %>%
  mutate(status = case_when(
    birds>0 & mammals==0 ~ "birds_only",
    birds==0 & mammals>0 ~ "mammals_only",
    birds>0 & mammals>0 ~ "birds_and_mammals",
    birds==0 & mammals==0 ~ "neither")) %>%
  with(.,table(status))




#############West Nile Virus
all_res_anno_wnv <- all_res_anno %>%
  filter(infect_type_virus_type == "west_nile_virus")

#Get sig counts for each gene for both birds and mammals
ens_sig_birds_mammals_wnv <- all_res_anno_wnv %>%
  filter(!is.na(ensembl_id_hs)) %>%
  mutate(sig_notsig = if_else(sig != 0, 1,0),up_reg = if_else(sig==1,1,0),down_reg = if_else(sig ==-1,1,0)) %>%
  group_by(ensembl_id_hs,clade) %>%
  summarize(n_sig = sum(sig_notsig)) %>%
  spread(clade,n_sig)

#Filter out missing genes in birds, make a table with counts in each category to get a sense of the data
ens_sig_birds_mammals_wnv %>%
  filter(!is.na(birds), !is.na(mammals)) %>%
  mutate(status = case_when(
    birds>0 & mammals==0 ~ "birds_only",
    birds==0 & mammals>0 ~ "mammals_only",
    birds>0 & mammals>0 ~ "birds_and_mammals",
    birds==0 & mammals==0 ~ "neither")) %>%
  with(.,table(status))





##########paramyxovirus #Prob not useful
all_res_anno_paramyxo<- all_res_anno %>%
  filter(infect_type_virus_type == "paramyxovirus")

#Get sig counts for each gene for both birds and mammals
ens_sig_birds_mammals_paramyxo <- all_res_anno_paramyxo %>%
  filter(!is.na(ensembl_id_hs)) %>%
  mutate(sig_notsig = if_else(sig != 0, 1,0),up_reg = if_else(sig==1,1,0),down_reg = if_else(sig ==-1,1,0)) %>%
  group_by(ensembl_id_hs,clade) %>%
  summarize(n_sig = sum(sig_notsig)) %>%
  spread(clade,n_sig)

#Filter out missing genes in birds, make a table with counts in each category to get a sense of the data
ens_sig_birds_mammals_paramyxo %>%
  filter(!is.na(birds), !is.na(mammals)) %>%
  mutate(status = case_when(
    birds>0 & mammals==0 ~ "birds_only",
    birds==0 & mammals>0 ~ "mammals_only",
    birds>0 & mammals>0 ~ "birds_and_mammals",
    birds==0 & mammals==0 ~ "neither")) %>%
  with(.,table(status))







############Ecoli
all_res_anno_ecoli<- all_res_anno %>%
  filter(infect_type_virus_type == "ecoli")

#Get sig counts for each gene for both birds and mammals
ens_sig_birds_mammals_ecoli <- all_res_anno_ecoli %>%
  filter(!is.na(ensembl_id_hs)) %>%
  mutate(sig_notsig = if_else(sig != 0, 1,0),up_reg = if_else(sig==1,1,0),down_reg = if_else(sig ==-1,1,0)) %>%
  group_by(ensembl_id_hs,clade) %>%
  summarize(n_sig = sum(sig_notsig)) %>%
  spread(clade,n_sig))

#Filter out missing genes in birds, make a table with counts in each category to get a sense of the data
ens_sig_birds_mammals_ecoli %>%
  filter(!is.na(birds), !is.na(mammals)) %>%
  mutate(status = case_when(
    birds>0 & mammals==0 ~ "birds_only",
    birds==0 & mammals>0 ~ "mammals_only",
    birds>0 & mammals>0 ~ "birds_and_mammals",
    birds==0 & mammals==0 ~ "neither")) %>%
  with(.,table(status))





