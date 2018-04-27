setwd("~/Dropbox/BirdImmuneGeneEvolution")
library(tidyverse)
library(myTAI)

#Read in enard et al. 2016 VIP and virus info
virus_names <- read_csv("07_input_virus_classification/enard_virus_names.csv")
vip_evi <- read_csv("07_input_virus_classification/enard_vip_evidence_only.csv",col_names = c("ensembl_hs","hgnc","evidence","host_species"))


#First get the taxonomy for each virus species used to identify VIPs, grab rank below kingdom, order, family, subfamily (if exists), and genus. Use NCBI taxonomy.
#Function to query the NCBI database, and return the fields of interest
test <- head(virus_names)

get_virus_tax <- function(virus_name,rank){
  tax <- taxonomy(virus_name)
  tax_rank <- tax[tax$rank == rank,"name"]
  if (length(tax_rank)>0){
    return(tax_rank)
  }
  else {
    return(NA)
  }
}

virus_names_anno <- virus_names %>%
  mutate(family = unlist(map(.$virus_name_corrected,function(x) get_virus_tax(x,"family"))))

virus_names_anno %>%
  distinct(family) %>%
  print(n=30)



