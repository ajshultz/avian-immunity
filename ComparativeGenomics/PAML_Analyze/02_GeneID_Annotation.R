setwd("~/Dropbox/BirdImmuneGeneEvolution")
library(tidyverse)
library(biomaRt)

#######################################################################################################################
#Assign NCBI numbers
#######################################################################################################################
#load previous results if run in previous session
load("01_output_processed_data/all_res_allgenes.RDat")

#Read in HOG IDs
hog_ids <- read.table("02_input_annotation_data/new_hog_list.txt")
colnames(hog_ids) <- c("HOG2_HogID","NCBI_ID","entrezgene","sp")

#Data cleanup, extract chicken gene IDs
hog_ids_gg <- hog_ids %>% tbl_df %>%
  separate(HOG2_HogID,sep = "_",into=c("drop","hog")) %>%
  mutate(entrezgene = as.character(entrezgene), sp = as.character(sp)) %>%
  filter(sp == "galGal") %>%
  dplyr::select(hog,entrezgene)

#Data cleanup, extract zebra finch gene IDs
hog_ids_zf <- hog_ids %>% tbl_df %>%
  separate(HOG2_HogID,sep = "_",into=c("drop","hog")) %>%
  mutate(entrezgene_zf = as.character(entrezgene), sp = as.character(sp)) %>%
  filter(sp == "taeGut") %>%
  dplyr::select(hog,entrezgene_zf)

#Merge Hog IDs and NCBI gene IDs (both chicken and zebra finch)
all_res_gene_ncbi <- all_res_gene %>%
  left_join(hog_ids_gg,by="hog") %>%
  left_join(hog_ids_zf,by="hog") %>%
  distinct(hog,.keep_all=TRUE)

all_res_sp_ncbi <- all_res_sp %>%
  left_join(hog_ids_gg,by="hog") %>%
  left_join(hog_ids_zf,by="hog") %>%
  distinct(hog,.keep_all=TRUE)

write_csv(all_res_gene_ncbi,path="02_output_annotated_data/raw_results_genetrees.csv")
write_csv(all_res_sp_ncbi,path="02_output_annotated_data/raw_results_sptrees.csv")

save(all_res_gene_ncbi,all_res_sp_ncbi,file="02_output_annotated_data/all_res_ncbi.Rdat")


####Add human NCBI gene IDs from biomaRt
ensembl = useMart("ensembl")

#Use correct datasets - sometimes throws an error for some reason saying the given dataset is not valid. However, it works if you try again once or twice
ensembl_gg = useDataset("ggallus_gene_ensembl",mart=ensembl)
ensembl_hs = useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl_zf = useDataset("tguttata_gene_ensembl",mart=ensembl)

#attributes = listAttributes(ensembl)
#filters = listFilters(ensembl)

#Get human ensembl IDs for homologs
human_ids_gg <- getBM(attributes=c("ensembl_gene_id","external_gene_name","hsapiens_homolog_ensembl_gene","hsapiens_homolog_orthology_confidence"),
                   filters=c("with_hsapiens_homolog"),
                   values=TRUE,
                   mart=ensembl_gg) %>% tbl_df

human_ids_zf <- getBM(attributes=c("ensembl_gene_id","external_gene_name","hsapiens_homolog_ensembl_gene","hsapiens_homolog_orthology_confidence"),
                    filters=c("with_hsapiens_homolog"),
                    values=TRUE,
                    mart=ensembl_zf) %>% tbl_df %>%
  dplyr::rename(ensembl_gene_id_zf = ensembl_gene_id)


#Get ensembl ID to NCBI entrezgene mapping for three different species
gg_ens_to_entrezgene <- getBM(attributes=c("ensembl_gene_id","entrezgene"),
  mart=ensembl_gg) %>%
  tbl_df %>%
  filter(!is.na(entrezgene)) %>%
  mutate(entrezgene = as.character(entrezgene))

zf_ens_to_entrezgene <- getBM(attributes=c("ensembl_gene_id","entrezgene"),
                              mart=ensembl_zf) %>%
  tbl_df %>%
  filter(!is.na(entrezgene)) %>%
  dplyr::rename(entrezgene_zf=entrezgene,ensembl_gene_id_zf=ensembl_gene_id) %>%
  mutate(entrezgene_zf = as.character(entrezgene_zf))

hs_ens_to_entrezgene <- getBM(attributes=c("ensembl_gene_id","entrezgene"),
  mart=ensembl_hs) %>%
  tbl_df %>%
  filter(!is.na(entrezgene)) %>%
  dplyr::rename(entrezgene_hs=entrezgene,hsapiens_homolog_ensembl_gene=ensembl_gene_id) %>%
  mutate(entrezgene_hs = as.character(entrezgene_hs))

#Add human IDs to chicken ensembl to entrezgene table

gg_trans_table <- gg_ens_to_entrezgene %>%
  left_join(human_ids_gg,by="ensembl_gene_id") %>%
  left_join(hs_ens_to_entrezgene,by="hsapiens_homolog_ensembl_gene")

zf_trans_table <- zf_ens_to_entrezgene %>%
  left_join(human_ids_zf,by="ensembl_gene_id_zf") %>%
  left_join(hs_ens_to_entrezgene,by="hsapiens_homolog_ensembl_gene")

#save biomart results, note saving disabled so that new annotations do not overwrite those used to conduct all analyses presented in paper.
#save(gg_ens_to_entrezgene, zf_ens_to_entrezgene, hs_ens_to_entrezgene, gg_trans_table, zf_trans_table, file="02_output_annotated_data/biomart_translation_tables.Rdat")

#######Load version from submitted paper, human annotations seemed to have changed a bit since then
load("02_output_annotated_data/biomart_translation_tables.Rdat")

##########Create combined dataset of all hog results, chicken and zebra finch annotations. Use the coalesce() function to replace missing chicken annotations with zebra finch. Note that this results in multiples in the case of several matches. For now keeping to make sure the highest possible number of matches with mammal dataset. Will get distinct values before analyses.
all_res_gene_zf_hs <- all_res_gene_ncbi %>%
  left_join(gg_trans_table,by="entrezgene") %>%
  left_join(zf_trans_table,by="entrezgene_zf") %>%
  mutate(entrezgene_hs = coalesce(entrezgene_hs.x,entrezgene_hs.y), ensembl_gene_id_hs = dplyr::coalesce(hsapiens_homolog_ensembl_gene.x,hsapiens_homolog_ensembl_gene.y))
write_delim(all_res_gene_zf_hs,path="02_output_annotated_data/raw_results_genetrees_all_annotations.csv",delim=",")

all_res_sp_zf_hs <- all_res_sp_ncbi %>%
  left_join(gg_trans_table,by="entrezgene") %>%
  left_join(zf_trans_table,by="entrezgene_zf") %>%
  mutate(entrezgene_hs = coalesce(entrezgene_hs.x,entrezgene_hs.y), ensembl_gene_id_hs = coalesce(hsapiens_homolog_ensembl_gene.x,hsapiens_homolog_ensembl_gene.y))
write_delim(all_res_sp_zf_hs,path="02_output_annotated_data/raw_results_sptrees_all_annotations.csv",delim=",")


save(all_res_gene_zf_hs,all_res_sp_zf_hs,file="02_output_annotated_data/all_res_zf_hs.Rdat")
