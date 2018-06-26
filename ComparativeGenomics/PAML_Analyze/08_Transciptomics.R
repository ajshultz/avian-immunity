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
#Note, have to create new bioproject IDs, as current does not include all conditions
all_res_birds <- all_res %>%
  separate(target_id,into=c("ensembl_id","drop"),fill = "right") %>%
  dplyr::select(-drop) %>%
  unite("bioproject_infect",c("bioproject","condition","infect"),sep = "-",remove=FALSE) %>%
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


#For all subsequent bioprojects, this function will create a clean dataset with only genes given as "signficiant" given a certain number of matches. User needs to specify which column contains the ensembl ids to group, and has the option to remove any treatments (specify vector with treatments (bioproject column) to remove). Also, if data are a time course with an additional overall comparison, user can specify whether to use the n_matches for all but an overall test column (e.g. time course with overall significance), which will require all other tests to be signficant > n_matches, as well as being signficant in the specified overall test column.
collapse_n_matches <- function(dataset,bioproject_number,n_matches,ensembl_id_column,remove_treatment=NULL,overall_test = "no",overall_column = NULL){
  
  #First remove treatment if necessary
  if (length(remove_treatment)>0 & !is.na(remove_treatment)) {
    dataset_filtered <- dataset %>%
      filter(!(bioproject %in% remove_treatment))
  } else {
    dataset_filtered <- dataset
  }
  
  if (overall_test == "no"){
    sig_values_fixed <- dataset_filtered %>%
      filter(bioproj_number == bioproject_number) %>%
      group_by(!!as.name(ensembl_id_column)) %>%
      summarize(beta_sig = sum(beta_sig), sig = sum(sig))  %>%
      mutate(singles_beta_sig = case_when(beta_sig >= n_matches ~ 1,
                                          beta_sig < n_matches & beta_sig >-n_matches ~ 0,
                                          beta_sig <= -n_matches ~ -1),
             singles_sig = case_when(sig >= n_matches ~ 1,
                                     sig < n_matches & sig >-n_matches ~ 0,
                                     sig <= -n_matches ~ -1)) %>%
      dplyr::select(-beta_sig,-sig) %>%
      ungroup()
  }
  if (overall_test == "yes") {
    dataset_nooverall <-  dataset_filtered %>%
      filter(bioproject != overall_column) %>%
      filter(bioproj_number == bioproject_number) %>%
      group_by(!!as.name(ensembl_id_column)) %>%
      summarize(beta_sig = sum(beta_sig), sig = sum(sig))  %>%
      mutate(singles_beta_sig = case_when(beta_sig >= n_matches ~ 1,
                                          beta_sig < n_matches & beta_sig >-n_matches ~ 0,
                                          beta_sig <= -n_matches ~ -1),
             singles_sig = case_when(sig >= n_matches ~ 1,
                                     sig < n_matches & sig >-n_matches ~ 0,
                                     sig <= -n_matches ~ -1)) %>%
      dplyr::select(-beta_sig,-sig) %>%
      ungroup()
    
    sig_values_fixed <- dataset_filtered %>%
      filter(bioproject == overall_column) %>%
      dplyr::select(!!as.name(ensembl_id_column),singles_beta_sig = beta_sig,singles_sig = sig) %>%
      bind_rows(dataset_nooverall) %>%
      group_by(!!as.name(ensembl_id_column)) %>%
      summarize(beta_sig = sum(singles_beta_sig), sig = sum(singles_sig))  %>%
      mutate(singles_beta_sig = case_when(beta_sig >= 2 ~ 1,
                                          beta_sig < 2 & beta_sig >-2 ~ 0,
                                          beta_sig <= -2 ~ -1),
             singles_sig = case_when(sig >= 2 ~ 1,
                                     sig < 2 & sig >-2 ~ 0,
                                     sig <= -2 ~ -1)) %>%
      dplyr::select(-beta_sig,-sig) %>%
      ungroup()
  }
  
  #Combine back with all_res_birds, remove duplicates and drop old columns
  all_res_clean <- dataset_filtered %>%
    filter(bioproj_number == bioproject_number) %>%
    left_join(sig_values_fixed,by=paste0(ensembl_id_column)) %>%
    dplyr::select(-beta_sig,-sig) %>%
    rename(beta_sig = singles_beta_sig,sig=singles_sig) %>%
    distinct(bioproj_number,!!as.name(ensembl_id_column),.keep_all=T)
  
  return(all_res_clean)
}

#Read in dtafraem with information on condition requirements to clean up data
bioprojects_cond_anno <- read_csv("08_inputs_transcriptomics/bioproject_cond_count_anno.csv")

bioprojects_cond_anno_birds <- bioprojects_cond_anno %>%
  filter(clade == "birds") %>%
  as.data.frame

#First, clean up bird-only datasets, so that all genes tested in our study can be compared to transcriptome enrichment. Iterate through each bioproject and create clean dataset wtih collapse_n_matches function, add to list and then bind together
bioproject_res_clean <- list()

for (i in 1:nrow(bioprojects_cond_anno_birds)){
  bioproject_res_clean[[i]] <- collapse_n_matches(all_res_birds_singles,bioproject_number=bioprojects_cond_anno_birds[i,"bioproj_number"], n_matches=bioprojects_cond_anno_birds[i,"n_required"], ensembl_id_column="ensembl_id_gg",remove_treatment=bioprojects_cond_anno_birds[i,"remove_treatment_column"],overall_test=bioprojects_cond_anno_birds[i,"overall_test"],overall_column=bioprojects_cond_anno_birds[i,"overall_column"])
}

all_res_birds_clean <- bind_rows(bioproject_res_clean)

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

#Get list of infectious agents to iterate through
infect_agents <- all_res_birds_clean %>%
  distinct(infect_type_virus_type) %>%
  pull

#Colors to use:
general_colors_plasmodium <- c("white","#88CCEE")
names(general_colors_plasmodium) <- c("not significant","significant")

general_colors_virus <- c("white","#44AA99")
names(general_colors_virus) <- c("not significant","significant")

general_colors_bacterium<- c("white","#332288")
names(general_colors_bacterium) <- c("not significant","significant")

color_list <- list(general_colors_plasmodium,general_colors_virus,general_colors_bacterium)
names(color_list) <- c("plasmodium","virus","bacterium")

color_vec <- as.tibble(color_list)[2,] %>% gather(agent,color) %>% pull(color)
names(color_vec) <- as.tibble(color_list)[2,] %>% gather(agent,color) %>% pull(agent)


#Create table to grab significance resutls
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

#Clean up actual results table:

infect_res_table_clean <- infect_res_table %>%
  as.tibble %>%
  mutate(estimate=round(as.numeric(estimate),digits=2),
         std.error=round(as.numeric(std.error),digits=2),
         z.value=round(as.numeric(z.value),digits=2),
         p.value=round(as.numeric(p.value),digits=4),
         prop_sig=round(as.numeric(prop_sig),digits=2))

write_csv(infect_res_table_clean,"08_output_transcriptomics/bird_transcriptome_sig_results_table.csv")

#########################################################################################################
################### Mammals Transcriptome dataset cleanup #############################################
#########################################################################################################
#Now, we will clean up bird and mammal data together - have to redo birds with human gene IDs

#First, need to get rid of extra info on end of target_ids, and create new column of sig given beta > 1 or -1
all_res_mammals <- all_res %>%
  separate(target_id,into=c("ensembl_id","drop"),fill = "right") %>%
  dplyr::select(-drop) %>%
  unite("bioproject_infect",c("bioproject","condition","infect"),sep = "-",remove=FALSE) %>%
  dplyr::select(-bioproject) %>%
  dplyr::rename(bioproject = bioproject_infect) %>%
  left_join(trans_table,by=c("ensembl_id" = "ensembl_gene_id")) %>%
  mutate(ensembl_id_hs = if_else(species=="HOM",ensembl_id,hsapiens_homolog_ensembl_gene)) %>%
  mutate(beta_TF = if_else(abs(b)>=1,1,0), beta_sig = sig*beta_TF) %>%
  separate(bioproject, into=c("bioproj_number"), extra="drop",remove = F,sep = "-") %>%
  left_join(infect_info,by=c("infect","species"))

#Now, want to make sure there are not duplicate gene IDs in each bioproject. For those that are, take the sum of beta_sig and sig, and change any >1 to 1 and <-1 to -1.
all_res_mammals_nodups <- all_res_mammals %>%
  group_by(bioproject,ensembl_id_hs) %>%
  summarize(beta_sig = sum(beta_sig), sig = sum(sig)) %>%
  mutate(singles_beta_sig = case_when(beta_sig >= 1 ~ 1,
                                      beta_sig == 0 ~ 0,
                                      beta_sig <= -1 ~ -1),
         singles_sig = case_when(sig >= 1 ~ 1,
                                 sig == 0 ~ 0,
                                 sig <= -1 ~ -1)) %>%
  dplyr::select(-beta_sig,-sig) %>%
  ungroup()

#Combine back with all_res_mammals, remove duplicates and drop old columns
all_res_mammals_singles <- all_res_mammals %>%
  left_join(all_res_mammals_nodups,by=c("bioproject","ensembl_id_hs")) %>%
  dplyr::select(-beta_sig,-sig) %>%
  rename(beta_sig = singles_beta_sig,sig=singles_sig) %>%
  distinct(bioproject,ensembl_id_hs,.keep_all=T)

#Now each specific dataset is cleaned up, but we need to clean up bioprojects that had more than 1 condition (timepoints, etc.)

#First, write a csv of how many conditions are in each bioproject (note that this was the base used to create the input for bioproject_cond_count_anno)
bioprojects_n_conditions_mammals <- all_res_mammals_singles %>%
  group_by(bioproj_number,bioproject) %>%
  distinct(bioproject) %>%
  group_by(bioproj_number) %>%
  summarize(n_cond = n()) %>%
  left_join(bioproject_info %>% distinct(bioproj_number,.keep_all=TRUE),by="bioproj_number") %>%
  write_csv("08_output_transcriptomics/bioproject_cond_count.csv")

#Write some basic stats for each bioproject + conditions
#Create a tibble with useful info for each bioproject
bioproject_info <- all_res_mammals %>%
  distinct(bioproject,.keep_all=T) %>%
  dplyr::select(bioproject,bioproj_number,infect,infect_type_specific,infect_type_virus_type,infect_type_general,clade)

all_res_mammals %>%
  filter(!is.na(sig)) %>%
  group_by(bioproject) %>%
  summarize(n_up = sum(sig==1),n_up_beta = sum(beta_sig==1), n_down = sum(sig==-1), n_down_beta = sum(beta_sig==-1), n_notsig = sum(sig==0), n_notsig_beta = sum(beta_sig ==0),n_ensID_hs = sum(!is.na(ensembl_id_hs))) %>%
  separate(bioproject,into=c("bioproj_number"),sep = "-",remove=F,extra = "drop") %>%
  left_join(bioproject_info %>% distinct(bioproj_number,.keep_all=TRUE), by="bioproj_number") %>%
  write_csv("08_output_transcriptomics/bioproject_basic_info.csv")

#Now clean up each bioproject according to criteria in bioproject_cond_anno
bioproject_res_clean_mammals <- list()
bioprojects_cond_anno_df <- bioprojects_cond_anno %>% as.data.frame

for (i in 1:nrow(bioprojects_cond_anno)){
  bioproject_res_clean_mammals[[i]] <- collapse_n_matches(all_res_mammals_singles,bioproject_number=bioprojects_cond_anno_df[i,"bioproj_number"], n_matches=bioprojects_cond_anno_df[i,"n_required"], ensembl_id_column="ensembl_id_hs",remove_treatment=bioprojects_cond_anno_df[i,"remove_treatment_column"],overall_test=bioprojects_cond_anno_df[i,"overall_test"],overall_column=bioprojects_cond_anno_df[i,"overall_column"])
}

all_res_mammals_clean <- bind_rows(bioproject_res_clean_mammals)

save(all_res_mammals_clean,file = "08_output_transcriptomics/mammals_birds_clean_transcriptomic_results.Rdat")


#Read in BUSTED significance results for bird and mammal combined dataset
imm <- read_csv("05_output_bird_mammal_comparison_results/bird_mammal_combined_dataset.csv")

#Combine transcriptome and selection results, clean up only relevant columns, remove any genes without mammal_q, bird_q,or beta_sig
trans_imm <- all_res_mammals_clean %>%
  left_join(imm,by=c("ensembl_id_hs" = "ensembl_gene_id_hs")) %>%
  dplyr::select(ensembl_id_hs,entrezgene_hs,bioproj_number,infect,species,infect_type_specific:sig,hog,sig_all,mammal_logp,bird_logp,mammal_q,bird_q) %>%
  filter(!is.na(mammal_q), !is.na(bird_q), !is.na(beta_sig)) 

#For getting bird and mammal up-reg and down-reg sig results, use qval<0.05
qval <- 0.05
#Specify which genes are sig in birds, mammals, or birds and mammals given a qvalue
trans_imm<- trans_imm %>%
  mutate(sig_birds = if_else(bird_q<=qval,1,0),
         sig_mammals = if_else(mammal_q<=qval,1,0),
         sig_both = if_else(bird_q<=qval & mammal_q<=qval,1,0))

#How many genes in each bioproject with significance results?
trans_imm %>%
  group_by(clade,infect_type_specific,bioproj_number) %>%
  summarize(n_genes = n(),n_up = sum(beta_sig==1), n_down=sum(beta_sig==-1), nsig_birds=sum(sig_birds), nsig_mammals=sum(sig_mammals), nsig_both=sum(sig_both)) %>%
  write_csv("08_output_transcriptomics/bioproj_counts_combo_sig_results.csv")

#Create tibble of only sig results, for later use
sig_results_simple <- trans_imm %>%
  distinct(ensembl_id_hs,.keep_all=TRUE) %>%
  dplyr::select(ensembl_id_hs,sig_birds,sig_mammals,sig_both)

#####################################################################################################
#Redo bird sig testing to see if results hold with Busted results and smaller set of genes
#####################################################################################################
#Create table to grab significance resutls
infect_res_table <- matrix(nrow=length(infect_agents)*3,ncol=8)
colnames(infect_res_table) <- c("infect_agent", "trans_response", "n", "estimate","std.error","z.value","p.value", "prop_sig")

infect_plots <- list()

start <- 1

for (i in 1:length(infect_agents)){
  infect_res_table[(start:(start+2)),1] <- infect_agents[i]
  infect_res <- trans_imm %>%
    filter(clade == "birds") %>%
    filter(infect_type_virus_type==infect_agents[i]) %>%
    group_by(ensembl_id_hs) %>%
    summarize(beta_sig_overall = sum(beta_sig), sig_overall = sum(sig)) %>%
    mutate(beta_sig = case_when(beta_sig_overall >= 1 ~ 1,
                                beta_sig_overall == 0 ~ 0,
                                beta_sig_overall <= -1 ~ -1),
           sig = case_when(sig_overall >= 1 ~ 1,
                           sig_overall == 0 ~ 0,
                           sig_overall <= -1 ~ -1)) %>%
    filter(ensembl_id_hs != "", !is.na(beta_sig)) %>%
    mutate(expr_sig = if_else(beta_sig !=0,1,0), up_reg = if_else(beta_sig==1,1,0), down_reg = if_else(beta_sig==-1,1,0), not_reg = if_else(beta_sig == 0,1,0)) %>%
    left_join(sig_results_simple)
  
  infect_res_table[start,2] <- "down"   
  try(infect_res_downreg_test <- infect_res %>%
        filter(beta_sig != 1) %>%
        glm(sig_birds ~ down_reg, family="binomial",data=.))
  try(infect_res_table[start,4:7] <- summary(infect_res_downreg_test)$coefficient[2,])
  
  infect_res_table[start+1,2] <- "none" 
  infect_res_table[start+1,4:7] <- NA
  
  infect_res_table[start+2,2] <- "up"   
  try(infect_res_upreg_test <- infect_res %>%
        filter(beta_sig != -1) %>%
        glm(sig_birds ~ up_reg, family="binomial",data=.))
  try(infect_res_table[start+2,4:7] <- summary(infect_res_upreg_test)$coefficient[2,])
  
  general_type <- bioproj_infect_info %>%
    filter(infect_type_virus_type == infect_agents[i]) %>%
    distinct(infect_type_general) %>%
    pull
  
  #Calculate proportion significant
  infect_res_table[start:(start+2),8] <- infect_res %>%
    with(.,table(sig_birds,beta_sig)) %>%
    as.tibble %>%
    mutate(beta_sig = case_when(beta_sig == -1 ~ "down",
                                beta_sig == 0 ~ "none",
                                beta_sig == 1 ~ "up")) %>%
    group_by(beta_sig) %>%
    mutate(n_cat = sum(n)) %>%
    spread(sig_birds,n) %>%
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
    ylim(0,1.1) +
    geom_text(aes(trans_response,label=n),nudge_y=0.08,size=4) +
    geom_text(aes(trans_response,label=sig_char),nudge_y=0.01,size=4)
  
  start <- start + 3 
}

names(infect_plots) <- infect_agents

birds_legend <- trans_imm %>%
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
ggsave("08_output_transcriptomics/bird_transcriptome_sig_figure_busted_q0.05_mammal_genes.pdf",height=7,width=11)

#Write table of results
infect_res_table_clean <- infect_res_table %>%
  as.tibble %>%
  mutate(estimate=round(as.numeric(estimate),digits=2),
         std.error=round(as.numeric(std.error),digits=2),
         z.value=round(as.numeric(z.value),digits=2),
         p.value=round(as.numeric(p.value),digits=4),
         prop_sig=round(as.numeric(prop_sig),digits=2))

write_csv(infect_res_table_clean,"08_output_transcriptomics/bird_transcriptome_sig_results_table_mammal_bird_gene_set.csv")



###########################################################################################################
#Do mammals show the same results?
#########################################################################################################
#Get list of infectious agents to iterate through
infect_agents_mammals <- all_res_mammals_clean %>%
  filter(clade == "mammals") %>%
  distinct(infect_type_virus_type) %>%
  pull

infect_res_table <- matrix(nrow=length(infect_agents_mammals)*3,ncol=8)
colnames(infect_res_table) <- c("infect_agent", "trans_response", "n", "estimate","std.error","z.value","p.value", "prop_sig")

infect_plots <- list()

start <- 1

for (i in 1:length(infect_agents_mammals)){
  infect_res_table[(start:(start+2)),1] <- infect_agents_mammals[i]
  infect_res <- trans_imm %>%
    filter(clade == "mammals") %>%
    filter(infect_type_virus_type==infect_agents_mammals[i]) %>%
    group_by(ensembl_id_hs) %>%
    summarize(beta_sig_overall = sum(beta_sig), sig_overall = sum(sig)) %>%
    mutate(beta_sig = case_when(beta_sig_overall >= 1 ~ 1,
                                beta_sig_overall == 0 ~ 0,
                                beta_sig_overall <= -1 ~ -1),
           sig = case_when(sig_overall >= 1 ~ 1,
                           sig_overall == 0 ~ 0,
                           sig_overall <= -1 ~ -1)) %>%
    filter(ensembl_id_hs != "", !is.na(beta_sig)) %>%
    mutate(expr_sig = if_else(beta_sig !=0,1,0), up_reg = if_else(beta_sig==1,1,0), down_reg = if_else(beta_sig==-1,1,0), not_reg = if_else(beta_sig == 0,1,0)) %>%
    left_join(sig_results_simple)
  
  infect_res_table[start,2] <- "down"   
  try(infect_res_downreg_test <- infect_res %>%
        filter(beta_sig != 1) %>%
        glm(sig_mammals ~ down_reg, family="binomial",data=.))
  try(infect_res_table[start,4:7] <- summary(infect_res_downreg_test)$coefficient[2,])
  
  infect_res_table[start+1,2] <- "none" 
  infect_res_table[start+1,4:7] <- NA
  
  infect_res_table[start+2,2] <- "up"   
  try(infect_res_upreg_test <- infect_res %>%
        filter(beta_sig != -1) %>%
        glm(sig_mammals ~ up_reg, family="binomial",data=.))
  try(infect_res_table[start+2,4:7] <- summary(infect_res_upreg_test)$coefficient[2,])
  
  general_type <- bioproj_infect_info %>%
    filter(infect_type_virus_type == infect_agents_mammals[i]) %>%
    distinct(infect_type_general) %>%
    pull
  
  #Calculate proportion significant
  infect_res_table[start:(start+2),8] <- infect_res %>%
    with(.,table(sig_mammals,beta_sig)) %>%
    as.tibble %>%
    mutate(beta_sig = case_when(beta_sig == -1 ~ "down",
                                beta_sig == 0 ~ "none",
                                beta_sig == 1 ~ "up")) %>%
    group_by(beta_sig) %>%
    mutate(n_cat = sum(n)) %>%
    spread(sig_mammals,n) %>%
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
    ggtitle(infect_agents_mammals[i]) +
    ylab("proportion significant") +
    xlab("transcriptional response") +
    ylim(0,1.1) +
    geom_text(aes(trans_response,label=n),nudge_y=0.08,size=4) +
    geom_text(aes(trans_response,label=sig_char),nudge_y=0.01,size=4)
  
  start <- start + 3 
}

names(infect_plots) <- infect_agents_mammals

mammals_legend <- trans_imm %>%
  ggplot(aes(infect_type_general,fill=factor(infect_type_general,levels=c("plasmodium","virus","bacterium")))) +
  geom_bar(position="fill") +
  scale_fill_manual(values = color_vec,name="pathogen type") +
  ylab("proportion") +
  ylim(1,2) +
  guides(fill=guide_legend(override.aes = list(fill=color_vec))) +
  theme(line = element_blank(),
        axis.text = element_blank(),
        title = element_blank())

plot_grid(infect_plots[["influenza"]],infect_plots[["coronavirus"]],infect_plots[["west_nile_virus"]],infect_plots[["ecoli"]],infect_plots[["mycoplasma"]],infect_plots[["plasmodium"]],mammals_legend,ncol=4)
ggsave("08_output_transcriptomics/mammal_transcriptome_sig_figure_busted_q0.05_mammal_genes.pdf",height=7,width=11)

#########################################################################################################
##################### Are genes up or downregulated in both birds and mammals more likely to be under selection in both clades?
#########################################################################################################
######First, is there more overlap than expected in which genes are up or down regulated in birds and mammals?

#How many genes in each bioproject with significance results?
trans_imm %>%
  group_by(clade,infect_type_specific,bioproj_number) %>%
  summarize(n_genes = n(),n_up = sum(beta_sig==1), n_down=sum(beta_sig==-1), nsig_birds=sum(sig_birds), nsig_mammals=sum(sig_mammals), nsig_both=sum(sig_both)) %>%
  write_csv("08_output_transcriptomics/bioproj_counts_combo_sig_q0.01_results.csv")


#Infectious agents tested in both clades
infect_agents_both <- c("influenza","west_nile_virus", "ecoli", "mycoplasma","plasmodium")

trans_res_table <- matrix(nrow=length(infect_agents_both)*2,ncol=8)
colnames(trans_res_table) <- c("infect_agent", "trans_response", "n", "n.expected", "p.value", "conf.int1", "conf.int2", "estimate.odds.ratio")

trans_res_list <- list()

trans_start <- 1

for (i in 1:length(infect_agents_both)){
  trans_res_table[(trans_start:(trans_start+1)),1] <- infect_agents_both[i]
  
  #Combine results across studies, then by clade, calculate overlaps
  infect_res <- trans_imm %>%
    filter(infect_type_virus_type==infect_agents_both[i]) %>%
    group_by(clade,ensembl_id_hs) %>%
    summarize(beta_sig_overall = sum(beta_sig), sig_overall = sum(sig)) %>%
    mutate(beta_sig = case_when(beta_sig_overall >= 1 ~ 1,
                                beta_sig_overall == 0 ~ 0,
                                beta_sig_overall <= -1 ~ -1),
           sig = case_when(sig_overall >= 1 ~ 1,
                           sig_overall == 0 ~ 0,
                           sig_overall <= -1 ~ -1)) %>%
    filter(ensembl_id_hs != "", !is.na(beta_sig)) %>%
    mutate(expr_sig = if_else(beta_sig !=0,1,0), up_reg = if_else(beta_sig==1,1,0), down_reg = if_else(beta_sig==-1,1,0), not_reg = if_else(beta_sig == 0,1,0)) %>%
    dplyr::select(clade,ensembl_id_hs,up_reg,down_reg,not_reg) 
  
  infect_res_birds <- infect_res %>%
    ungroup() %>%
    filter(clade=="birds") %>%
    dplyr::select(ensembl_id_hs,up_reg_birds=up_reg,down_reg_birds=down_reg,not_reg_birds=not_reg)
  
  infect_res_mammals <- infect_res %>%
    ungroup() %>%
    filter(clade=="mammals") %>%
    dplyr::select(ensembl_id_hs,up_reg_mammals=up_reg,down_reg_mammals=down_reg,not_reg_mamamls=not_reg)

  #Get joint transcriptome response vars from birds and mammals, combine with significance results
  infect_res <- infect_res %>%
    ungroup() %>%
    group_by(ensembl_id_hs) %>%
    summarize(up_reg_both = sum(up_reg), down_reg_both = sum(down_reg), not_reg_both = sum(not_reg)) %>%
    mutate(up_reg_both = if_else(up_reg_both == 2,1,0), down_reg_both = if_else(down_reg_both==2,1,0), not_reg_both = if_else(not_reg_both==2,1,0)) %>%
    left_join(infect_res_birds) %>%
    left_join(infect_res_mammals) %>%
    filter(!is.na(up_reg_birds), !is.na(up_reg_mammals)) %>%
    mutate(infect_agent = infect_agents_both[i])
  
  trans_res_list[[i]] <- infect_res
  
  #Is there a signficant overlap of genes up-regulated in birds and mammals?
  trans_res_table[trans_start,2] <- "up"
  trans_res_table[trans_start,3] <- infect_res %>%
    filter(down_reg_both != 1) %>%
    with(.,table(up_reg_mammals,up_reg_birds)) %>%
    .[2,2]
  trans_res_table[trans_start,5:8] <- infect_res %>%
    filter(down_reg_both != 1) %>%
    with(.,table(up_reg_mammals,up_reg_birds)) %>%
    fisher.test %>%
    unlist %>%
    .[1:4]
  
  #What is the expected number of genes up regulated in both?
  trans_res_table[trans_start,4] <- infect_res %>%
    filter(down_reg_both != 1) %>%
    with(.,table(up_reg_mammals,up_reg_birds)) %>%
    chisq.test() %>%
    .$expected %>%
    unlist %>%
    .[2,2]
  
  #Is there a signficant overlap of genes down-regulated in birds and mammals?
  trans_res_table[trans_start+1,2] <- "down"
  trans_res_table[trans_start+1,3] <- infect_res %>%
    filter(up_reg_both != 1) %>%
    with(.,table(down_reg_mammals,down_reg_birds)) %>%
    .[2,2]
  trans_res_table[trans_start+1,5:8] <- infect_res %>%
    filter(up_reg_both != 1) %>%
    with(.,table(down_reg_mammals,down_reg_birds)) %>%
    fisher.test %>%
    unlist %>%
    .[1:4]
  
  #What is the expected number of genes down regulated in both?
  trans_res_table[trans_start+1,4] <- infect_res %>%
    filter(up_reg_both != 1) %>%
    with(.,table(down_reg_mammals,down_reg_birds)) %>%
    chisq.test() %>%
    .$expected %>%
    unlist %>%
    .[2,2]
  
  trans_start = trans_start + 2
}

#Clean up results and save  
trans_res_table_clean <- trans_res_table %>%
  as.tibble %>%
  mutate(n_both = as.numeric(n), n_both_expected = round(as.numeric(n.expected),digits=1), p.value = round(as.numeric(p.value),digits=5), odds_ratio = round(as.numeric(estimate.odds.ratio),digits=2), upper_bound = round(as.numeric(conf.int2),digits=2), lower_bound = round(as.numeric(conf.int1),digits=2)) %>%
  dplyr::select(infect_agent,trans_response,n_both,n_both_expected,p.value,odds_ratio,lower_bound,upper_bound)
  
trans_res_table_clean %>%
  write_csv("08_output_transcriptomics/birds_mammals_transcriptome_odds_ratio_table.csv")

#Odds ratio plot
trans_res_table_clean %>%
  filter(infect_agent != "west_nile_virus") %>%
  ggplot(aes(infect_agent,odds_ratio,col=trans_response)) +
  geom_point(size=4,position=position_dodge(width=0.9)) +
  geom_linerange(aes(infect_agent,ymin=lower_bound,ymax=upper_bound,col=trans_response),size=2,position=position_dodge(width=0.9)) +
  geom_hline(aes(yintercept = 1),size=2,linetype="dashed",col="black") +
  xlab("infectious agent") +
  ylab("odds ratio") +
  ylim(0,10) +
  scale_color_manual(values=c("up"="#CC6677","down"="#332288"),name="response")
ggsave("08_output_transcriptomics/birds_mammals_transcriptome_odds_ratio_plot.pdf",width=7,height=4)
  
#Combine transcriptome results into a single tibble, can group or split by infect_agent column
trans_res_all <- bind_rows(trans_res_list)


################
#Now to add in significance results, identify genes that are selected in both, and up regulated or down regulated in both. Will test different q-values
#Recreating the trans_imm data structure to allow for different q values
qvals <- c(0.1,0.01,0.001,0.0001)

sig_res_both_list <- list()

for (i in 1:length(qvals)){
  #Specify which genes are sig in birds, mammals, or birds and mammals given a qvalue
  trans_imm<- trans_imm %>%
  mutate(sig_birds = if_else(bird_q<=qvals[i],1,0),
          sig_mammals = if_else(mammal_q<=qvals[i],1,0),
          sig_both = if_else(bird_q<=qval & mammal_q<=qvals[i],1,0))
  
  #Create simplified list to join to transcriptome results
  sig_results_simple <- trans_imm %>%
    distinct(ensembl_id_hs,.keep_all=TRUE) %>%
    dplyr::select(ensembl_id_hs,sig_birds,sig_mammals,sig_both)
  
  #Combine transcriptome and significance results
  trans_res_sig <- trans_res_all %>%
    left_join(sig_results_simple,by="ensembl_id_hs") 
  
  sig_res_both_list[[i]] <- trans_res_sig %>%
    filter(sig_both == 1, up_reg_both == 1 | down_reg_both == 1) %>%
    mutate(qval=qvals[i])
}

sig_res_both_all <- bind_rows(sig_res_both_list) %>%
  dplyr::select(qval,infect_agent,ensembl_id_hs,up_reg_both,down_reg_both,sig_both)

#Add gene names, write to file
sig_res_both_all %>%
  left_join(gg_trans_table,by=c("ensembl_id_hs"="hsapiens_homolog_ensembl_gene")) %>%
  dplyr::select(qval:sig_both,external_gene_name,entrezgene,entrezgene_hs) %>%
  write_csv("08_output_transcriptomics/birds_mammals_sig_both_reg_both_genes.csv")


####################
#Are genes up or down regulated in birds, more likley to be under selection in both lineages?
infect_res_table <- matrix(nrow=length(infect_agents_both)*6,ncol=8)
colnames(infect_res_table) <- c("infect_agent", "trans_response", "n", "test","estimate","std.error","z.value","p.value")

#Test q<0.05
qval <- 0.05

trans_imm<- trans_imm %>%
  mutate(sig_birds = if_else(bird_q<=qval,1,0),
         sig_mammals = if_else(mammal_q<=qval,1,0),
         sig_both = if_else(bird_q<=qval & mammal_q<=qval,1,0))

#Create simplified list to join to transcriptome results
sig_results_simple <- trans_imm %>%
  distinct(ensembl_id_hs,.keep_all=TRUE) %>%
  dplyr::select(ensembl_id_hs,sig_birds,sig_mammals,sig_both)
  
start <- 1

for (i in 1:length(infect_agents_both)){
  
  infect_res_table[start:(start+5),1] <- infect_agents_both[i]
  infect_res_table[start:(start+5),3] <- nrow(trans_res_list[[i]])
  
  infect_res_table[start:(start+2),2] <- "down"   
  try(infect_res_downreg_test <- trans_res_list[[i]] %>%
        left_join(sig_results_simple) %>%
        filter(up_reg_both != 1) %>%
        glm(sig_birds ~ sig_mammals*down_reg_both, family="binomial",data=.))
  try(infect_res_table[start:(start+1),4] <- rownames(summary(infect_res_downreg_test)$coefficient[2:3,]))
  try(infect_res_table[start:(start+1),5:8] <- summary(infect_res_downreg_test)$coefficient[2:3,])
  try(infect_res_table[start:(start+2),4] <- rownames(summary(infect_res_downreg_test)$coefficient[2:4,]))
  try(infect_res_table[start:(start+2),5:8] <- summary(infect_res_downreg_test)$coefficient[2:4,])

  infect_res_table[(start+3):(start+5),2] <- "up"   
  try(infect_res_upreg_test <- trans_res_list[[i]] %>%
        left_join(sig_results_simple) %>%
        filter(down_reg_birds != 1) %>%
        glm(sig_birds ~ sig_mammals*up_reg_birds, family="binomial",data=.))
  try(infect_res_table[(start+3):(start+4),4] <- rownames(summary(infect_res_upreg_test)$coefficient[2:3,]))
  try(infect_res_table[(start+3):(start+4),5:8] <- summary(infect_res_upreg_test)$coefficient[2:3,])
  try(infect_res_table[(start+3):(start+5),4] <- rownames(summary(infect_res_upreg_test)$coefficient[2:4,]))
  try(infect_res_table[(start+3):(start+5),5:8] <- summary(infect_res_upreg_test)$coefficient[2:4,])
  
  start <- start + 6
  }

#Clean up and write to file
infect_res_table %>%
  as.data.frame %>%
  as.tibble %>%
  mutate(estimate = round(as.numeric(as.character(estimate)),digits=2), std.error = round(as.numeric(as.character(std.error)),digits=2), z.value = round(as.numeric(as.character(z.value)),digits=2), p.value = round(as.numeric(as.character(p.value)),digits=5)) %>%
  write_csv("08_output_transcriptomics/birds_sig_trans_resp_and_sig_mammals_tests_q0.05.csv")
  


#################
#Now we want to compare the actual beta values between bird and mammal datasets.

#First, we need to standardize the beta values with z-normalization
trans_imm <- all_res_mammals_clean %>%
  group_by(bioproj_number) %>%
  mutate(b_scaled = (b-mean(b,na.rm=TRUE))/sd(b,na.rm=TRUE)) %>%
  ungroup() %>%
  left_join(imm,by=c("ensembl_id_hs" = "ensembl_gene_id_hs")) %>%
  dplyr::select(ensembl_id_hs,entrezgene_hs,bioproj_number,infect,species,infect_type_specific:sig,b_scaled,hog,sig_all,mammal_logp,bird_logp,mammal_q,bird_q) %>%
  filter(!is.na(mammal_q), !is.na(bird_q), !is.na(beta_sig)) 

#For getting bird and mammal up-reg and down-reg sig results, use qval<0.01
qval <- 0.05
#Specify which genes are sig in birds, mammals, or birds and mammals given a qvalue
trans_imm<- trans_imm %>%
  mutate(sig_birds = if_else(bird_q<=qval,1,0),
         sig_mammals = if_else(mammal_q<=qval,1,0),
         sig_both = if_else(bird_q<=qval & mammal_q<=qval,1,0))

#Create simplified list to join to transcriptome results, adding in sig_none and bird_sig_only
sig_results_simple <- trans_imm %>%
  distinct(ensembl_id_hs,.keep_all=TRUE) %>%
  dplyr::select(ensembl_id_hs,sig_birds,sig_mammals,sig_both) %>%
  mutate(sig_none = if_else(sig_birds == 0 & sig_mammals == 0, 1, 0),
         sig_birds_only = if_else(sig_birds ==1 & sig_mammals == 0,1,0),
         sig_cat = case_when(sig_birds == 0 & sig_mammals == 0 ~ "not_sig",
                             sig_birds == 0 & sig_mammals == 1 ~ "not_sig",
                             sig_birds == 1 & sig_mammals == 0 ~ "birds_only",
                             sig_birds == 1 & sig_mammals ==1 ~ "birds_and_mammals"))

#Now, for each clade and each infect_type_virus_type, take biggest (absolute value) beta score from replicate studies.
trans_infect_type <- trans_imm %>%
  group_by(clade,infect_type_virus_type,ensembl_id_hs) %>%
  summarize(max_b_scaled = max(abs(b_scaled)), b_dir = sum(b_scaled),beta_sig =  sum(beta_sig), sig = sum(sig)) %>%
  mutate(b_dir=case_when(b_dir>0~1,b_dir<0~-1,b_dir==0~0),
         beta_sig=case_when(beta_sig>0~1,beta_sig<0~-1,beta_sig==0~0),
         sig=case_when(sig>0~1,sig<0~-1,sig==0~0)) %>%
  ungroup()

harmMean <- function(x, removeNA=TRUE) {
  if (removeNA) {
    y=x[!is.na(x)]
    length(y)/sum(1/y)
  } else {
    length(x)/sum(1/x)
  }
}

sig_trans_infect_type <- trans_infect_type %>%
  group_by(infect_type_virus_type,ensembl_id_hs) %>%
  summarize(b_m_harm_mean_b = harmMean(max_b_scaled), b_m_b_dir_sum = sum(b_dir), b_m_beta_sig = sum(beta_sig), b_m_sig = sum(sig)) %>%
  filter(infect_type_virus_type %in% infect_agents_both) %>%
  ungroup() %>%
  left_join(sig_results_simple)


sig_trans_infect_type %>%
  #mutate(b_m_harm_mean_b=b_m_harm_mean_b*b_m_sig) %>%
  #filter(b_m_harm_mean_b!=0) %>%
  ggplot(aes(b_m_harm_mean_b,fill=factor(sig_birds))) +
  geom_density(alpha=0.5) +
  facet_wrap(~infect_type_virus_type)

sig_trans_infect_type %>%
  #mutate(b_m_harm_mean_b=b_m_harm_mean_b*b_m_sig) %>%
  #filter(b_m_harm_mean_b!=0) %>%
  ggplot(aes(b_m_harm_mean_b,fill=factor(sig_both))) +
  geom_density(alpha=0.5) +
  facet_wrap(~infect_type_virus_type)

sig_trans_infect_type %>%
  #mutate(b_m_harm_mean_b=b_m_harm_mean_b*b_m_sig) %>%
  #filter(b_m_harm_mean_b!=0) %>%
  ggplot(aes(b_m_harm_mean_b,fill=factor(sig_cat))) +
  geom_density(alpha=0.5) +
  facet_wrap(~infect_type_virus_type)

sig_trans_infect_type %>%
  #filter(b_m_sig!=0) %>%
  mutate(b_m_harm_mean_b = abs(b_m_harm_mean_b)) %>%
  ggplot(aes(sig_cat,b_m_harm_mean_b,fill=sig_cat)) +
  geom_boxplot(notch=TRUE,outlier.shape=NA,coef=1) +
  coord_cartesian(ylim=c(0,2.5)) +
  scale_fill_manual(guide=F,values=c("birds_and_mammals"="#CC6677", "birds_only"="#DDCC77", "not_sig"="#88CCEE")) +
  scale_x_discrete(labels=c("birds_and_mammals"="birds and mammals","birds_only"="birds only","not_sig"="not significant")) +
  theme(axis.text.x=element_text(angle=90))+
  ylab("abs(harmonic mean beta)")+
  xlab("significance category") +
  facet_grid(.~infect_type_virus_type) +
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),
        strip.background=element_rect(fill="white"))
ggsave("08_output_transcriptomics/birds_mammals_trans_sig_boxplot.pdf",width=10,height=6)

#Use a Mann-Whitney U-test to see if means are different between sig in birds and mammals and not sig
both_vs_not_trans_sig_res <- matrix(nrow=length(infect_agents_both),ncol=4)
birds_only_vs_not_trans_sig_res <- matrix(nrow=length(infect_agents_both),ncol=4)
both_vs_birds_only_trans_sig_res <- matrix(nrow=length(infect_agents_both),ncol=4)

for (i in 1:length(infect_agents_both)){
  both_vs_not_trans_sig_res[i,1] <- infect_agents_both[i]
  both_vs_not_trans_sig_res[i,2] <- "birds_and_mammals vs. not_sig"
  both_vs_not_trans_sig_res[i,3:4] <- sig_trans_infect_type %>%
    filter(infect_type_virus_type == infect_agents_both[i]) %>%
    #filter(b_m_sig != 0) %>%
    filter(sig_birds_only != 1) %>%
    spread(sig_both,b_m_harm_mean_b) %>%
    with(., wilcox.test(`0`,`1`)) %>%
    unlist() %>%
    .[c("statistic.W","p.value")]
  
  birds_only_vs_not_trans_sig_res[i,1] <- infect_agents_both[i]
  birds_only_vs_not_trans_sig_res[i,2] <- "birds_only vs. not_sig"
  birds_only_vs_not_trans_sig_res[i,3:4] <- sig_trans_infect_type %>%
    filter(infect_type_virus_type == infect_agents_both[i]) %>%
    #filter(b_m_sig != 0) %>%
    filter(sig_both != 1) %>%
    spread(sig_birds_only,b_m_harm_mean_b) %>%
    with(., wilcox.test(`0`,`1`)) %>%
    unlist() %>%
    .[c("statistic.W","p.value")]
  
  both_vs_birds_only_trans_sig_res[i,1] <- infect_agents_both[i]
  both_vs_birds_only_trans_sig_res[i,2] <- "birds_and_mammals vs. birds_only"
  both_vs_birds_only_trans_sig_res[i,3:4] <- sig_trans_infect_type %>%
    filter(infect_type_virus_type == infect_agents_both[i]) %>%
    #filter(b_m_sig != 0) %>%
    filter(sig_both != 0 | sig_birds_only != 0) %>%
    spread(sig_birds_only,b_m_harm_mean_b) %>%
    with(., wilcox.test(`0`,`1`)) %>%
    unlist() %>%
    .[c("statistic.W","p.value")]
}
colnames(both_vs_not_trans_sig_res) <- c("infect_agent","comparison","W","p.value")
colnames(birds_only_vs_not_trans_sig_res) <- c("infect_agent","comparison","W","p.value")
colnames(both_vs_birds_only_trans_sig_res) <- c("infect_agent","comparison","W","p.value")

both_vs_not_trans_sig_res <- both_vs_not_trans_sig_res %>%
  as.data.frame %>% as.tibble %>%
  mutate(W = as.numeric(as.character(W)), p.value = round(as.numeric(as.character(p.value)),digits=4))
birds_only_vs_not_trans_sig_res <- birds_only_vs_not_trans_sig_res %>%
  as.data.frame %>% as.tibble %>%
  mutate(W = as.numeric(as.character(W)), p.value = round(as.numeric(as.character(p.value)),digits=4))
both_vs_birds_only_trans_sig_res <- both_vs_birds_only_trans_sig_res %>%
  as.data.frame %>% as.tibble %>%
  mutate(W = as.numeric(as.character(W)), p.value = round(as.numeric(as.character(p.value)),digits=4))

all_trans_sig_res <- bind_rows(both_vs_not_trans_sig_res,birds_only_vs_not_trans_sig_res,both_vs_birds_only_trans_sig_res)  

all_trans_sig_res %>%
  write_csv("08_output_transcriptomics/birds_mammals_trans_sig_both_wilcox_test_results.csv")
