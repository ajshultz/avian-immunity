setwd("~/Dropbox/BirdImmuneGeneEvolution")
library(tidyverse)
library(myTAI)

#Read in enard et al. 2016 VIP and virus info
virus_names <- read_csv("07_input_virus_classification/enard_virus_names.csv")
vip_evi <- read_csv("07_input_virus_classification/enard_vip_evidence_only.csv",col_names = c("ensembl_hs","hgnc","evidence","host_species","comment"),skip = 1)
virus_trans_anno <- read_csv("07_input_virus_classification/viral_family_transmission_evidence.csv")


#First get the taxonomy for each virus species used to identify VIPs, grab rank below kingdom, order, family, subfamily (if exists), and genus. Use NCBI taxonomy.
#Function to query the NCBI database, and return the fields of interest

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


#Add protein evidence if no dna or rna
add_protein <- function(data){
  split_data <- unlist(strsplit(data,split = "-"))
  if (length(split_data) == 3) {
    if (split_data[2] != "dna" & split_data[2] != "rna" & split_data[2] != "vir") {
      new_split_data <- c(split_data[1],"protein",split_data[2:3])
      new_data <- paste(new_split_data,collapse = "-")
      return(new_data)
    }
    else{
      return(data)
    }
  }
  if (length(split_data) == 4){
    if (split_data[2] != "dna" & split_data[2] != "rna" & split_data[2] != "vir") {
      new_split_data <- c(split_data[1],"protein",split_data[2:4])
      new_data <- paste(new_split_data,collapse = "-")
      return(new_data) 
    }
    else{
      return(data)
    }
  }
  return(data)
}

#Separate out results, add protein evidence where missing, change two agents to separate lines when present
vip_evi_proc <- vip_evi %>%
  separate_rows(evidence,sep = ",") %>%
  mutate(new_evidence = unlist(map(.$evidence, add_protein))) %>%
  dplyr::select(-evidence) %>%
  separate(new_evidence,sep = "-",into=c("pubmed","molecule","agent1","agent2","type"),fill="right") %>%
  mutate(new_type = coalesce(type,agent2)) %>%
  mutate(agent2 = if_else(
    !(agent2 %in% c("dsDNA","dsDNART","ssRNART","ssRNA","ssRNA","ssDNA","dsRNA")),true=agent2,false=NA_character_)) %>%
  gather(ver,agent,agent1:agent2) %>%
  filter(!is.na(agent)) %>%
  dplyr::select(-ver,-type)

#Join evidence and virus names tibbles
vip_evi_proc_anno <- vip_evi_proc %>%
  left_join(virus_names_anno,by=c("agent"="abbreviation"))

#Write counts of genes associated with each family:
vip_evi_proc_anno %>%
  group_by(family) %>%
  summarize(count=n()) %>%
  write_csv("07_output_virus_classification/vip_family_counts.csv")

#Write counts of genes associated with each agent
vip_evi_proc_anno %>%
  group_by(agent) %>%
  summarize(count=n()) %>%
  write_csv("07_output_virus_classification/vip_agent_counts.csv")

#Clean up object, join to family transmission data
vip_evi_proc_anno <- vip_evi_proc_anno %>%
  dplyr::select(ensembl_hs,host_species,agent,virus_name = virus_name_corrected,family,type=new_type) %>%
  left_join(virus_trans_anno,by=c("family"))

#Save object
vip_evi_proc_anno %>%
  write_csv("07_output_virus_classification/vip_evidence_family_annotation.csv")

#Read in bird and mammal selection data
imm <- read_csv("05_output_bird_mammal_comparison_results/bird_mammal_combined_dataset.csv")

#Combine bird and mammal selection data and vip annotation
imm <- imm %>%
  left_join(vip_evi_proc_anno,by=c("ensembl_gene_id_hs" = "ensembl_hs"))

#How many genes are present in each viral family (remove family duplicates first)
imm %>%
  filter(!is.na(family)) %>%
  distinct(hog,family) %>%
  group_by(family) %>%
  summarize(count=n()) %>%
  print(n=30) %>%
  summarize(total_count=sum(count))

#How many genes are present in each transmission category  and viral type (remove duplicates first)
imm %>%
  filter(!is.na(trans_evidence)) %>%
  distinct(hog,trans_evidence) %>%
  group_by(trans_evidence) %>%
  summarize(count=n())

imm %>%
  filter(!is.na(type)) %>%
  distinct(hog,type) %>%
  group_by(type) %>%
  summarize(count=n())

#Make vector of families with at least 10 genes in the dataset
families_to_test <- imm %>%
  filter(!is.na(family)) %>%
  distinct(hog,family) %>%
  group_by(family) %>%
  summarize(count=n()) %>%
  filter(count>10) %>%
  pull(family)

#Make vector of transfer evidence type
trans_evidence_types <- imm %>%
  filter(!is.na(trans_evidence)) %>%
  distinct(trans_evidence) %>%
  pull()

#Make a vector of overall viral type
viral_type <-  imm %>%
  filter(!is.na(type)) %>%
  distinct(type) %>%
  pull()

imm %>%
  filter(!is.na(type)) %>%
  distinct(type,family) %>%
  print(n=30)

#Create column to identify non-immune genes (not vips, bips, or pips)
imm <- imm %>% mutate(non_immune = as.logical(vip == FALSE & bip == FALSE & pip == FALSE))

#######################################################################################################################
##################### First, look for enrichment in different classes of evidence types #####################
#######################################################################################################################


#q values to consider
qvals <- c(0.1,0.01,0.001,0.0001)

#List to capture results across q values
qval_res_trans_list <- list()

for (i in 1:length(qvals)){
  
  comp_propsig_trans <- matrix(nrow=(length(trans_evidence_types)+1),ncol=18)
  comp_propsig_trans[,1] <- c(trans_evidence_types,"all_genes")
  
  for (j in 1:length(trans_evidence_types)){
    #Create a vector of true and false for that evidence type for each hog, and remove duplicates
    trans_evi_type_hogs = imm %>%
      filter(trans_evidence == trans_evidence_types[j]) %>%
      pull(hog)
    imm_test <- imm %>%
      distinct(hog,.keep_all=TRUE)
      mutate(trans_type = if_else(hog %in% trans_evi_type_hogs,TRUE,FALSE))
    
    #Calculate fisher's exact test for prop selected in both birds and mammals, and get prop sig genes
    comp_propsig_trans[j,2:5] <- imm_test %>%
      with(.,table(trans_type,mammal_q < qvals[i] & bird_q < qvals[i])) %>% fisher.test %>% unlist %>% .[1:4]
    comp_propsig_trans[j,6] <- imm_test %>%
      filter(trans_type == TRUE) %>%
      with(.,prop.table(table(mammal_q < qvals[i] & bird_q < qvals[i]))) %>% .[2]
    
    #Calculate fisher's exact test for prop selected in mammals only, and get prop sig genes
    comp_propsig_trans[j,7:10] <- imm_test %>%
      with(.,table(trans_type,mammal_q < qvals[i])) %>% fisher.test %>% unlist %>% .[1:4]
    comp_propsig_trans[j,11] <- imm_test %>%
      filter(trans_type == TRUE) %>%
      with(.,prop.table(table(mammal_q < qvals[i]))) %>% .[2]
    
    #Calculate fisher's exact test for prop selected in birds only, and get prop sig genes
    comp_propsig_trans[j,12:15] <- imm_test %>%
      with(.,table(trans_type,bird_q < qvals[i])) %>% fisher.test %>% unlist %>% .[1:4]
    comp_propsig_trans[j,16] <- imm_test %>%
      filter(trans_type == TRUE) %>%
      with(.,prop.table(table(bird_q < qvals[i]))) %>% .[2]
    
    #Number of genes in each category:
    comp_propsig_trans[j,17] <- imm_test %>% filter(trans_type == TRUE) %>% summarize(n()) %>% pull
    }
  
  #Add prop sig all genes
  comp_propsig_trans[j+1,6] <- imm_test %>%
    with(.,prop.table(table(mammal_q < qvals[i] & bird_q < qvals[i]))) %>% .[2]
  comp_propsig_trans[j+1,11] <- imm_test %>%
    with(.,prop.table(table(mammal_q < qvals[i]))) %>% .[2]
  comp_propsig_trans[j+1,16] <- imm_test %>%
    with(.,prop.table(table(bird_q < qvals[i]))) %>% .[2]
  comp_propsig_trans[j+1,17] <- imm_test %>% summarize(n()) %>% pull
  
  comp_propsig_trans[,18] <- qvals[i]
  
  colnames(comp_propsig_trans) <- c("class","p.value_both","conf.int1_both","conf.int2_both","estimated.odds.ratio_both","prop_sel_both","p.value_mammals","conf.int1_mammals","conf.int2_mammals","estimated.odds.ratio_mammals","prop_sel_mammals","p.value_birds","conf.int1_birds","conf.int2_birds","estimated.odds.ratio_birds","prop_sel_birds","n.genes","qval")
  
  #Clean up, select relevant columns
  comp_propsig_trans_clean <- comp_propsig_trans %>%
    as.tibble %>%
    mutate(p.value_both=round(as.numeric(p.value_both),digits = 4),
          odds.ratio_both=round(as.double(estimated.odds.ratio_both),digits=2),
          prop.sel_both=round(as.numeric(prop_sel_both),digits=3),
          conf.lower_both=round(as.numeric(conf.int1_both),digits=3),
          conf.upper_both=round(as.numeric(conf.int2_both),digits=3),
          p.value_mammals=round(as.numeric(p.value_mammals),digits = 4),
          odds.ratio_mammals=round(as.numeric(estimated.odds.ratio_mammals),digits=2),
          prop.sel_mammals=round(as.numeric(prop_sel_mammals),digits=3),
          conf.lower_mammals=round(as.numeric(conf.int1_mammals),digits=3),
          conf.upper_mammals=round(as.numeric(conf.int2_mammals),digits=3),
          p.value_birds=round(as.numeric(p.value_birds),digits = 4),
          odds.ratio_birds=round(as.numeric(estimated.odds.ratio_birds),digits=2),
          prop.sel_birds=round(as.numeric(prop_sel_birds),digits=3),
          conf.lower_birds=round(as.numeric(conf.int1_birds),digits=3),
          conf.upper_birds=round(as.numeric(conf.int2_birds),digits=3)) %>%
    dplyr::select(qval,class,n.genes,odds.ratio_both,conf.lower_both,conf.upper_both,p.value_both,prop.sel_both,odds.ratio_mammals,conf.lower_mammals,conf.upper_mammals,p.value_mammals,prop.sel_mammals,odds.ratio_birds,conf.lower_birds,conf.upper_birds,p.value_birds,prop.sel_birds)
  
  
  qval_res_trans_list[[i]] <- comp_propsig_trans_clean
}

qval_res_trans <- qval_res_trans_list %>% bind_rows

write_csv(qval_res_trans,path="07_output_virus_classification/bird_mammals_evidence_type_enrichment_table.csv")


#Reshape into long format for plotting
qval_res_trans <- qval_res_trans %>%
  gather(stat,value,-qval,-class,-n.genes) %>%
  separate(stat,into = c("stat","comparison"),sep="_")

#Create vector of asterisks to include in plots:
qval_res_trans_forplotting <-   qval_res_trans %>%
  #filter(class != "all_genes") %>%
  filter(class != "only_mammals") %>%
  spread(stat,value) %>%
  mutate(sig = case_when(
    is.na(p.value) ~ "",
    p.value > 0.05 ~ "",
    p.value <= 0.05 & p.value > 0.01 ~ "*",
    p.value <= 0.01 & p.value > 0.001 ~ "**",
    p.value <= 0.001 ~ "***"
  ))


prop_sel_plot <- qval_res_trans_forplotting %>%
  #filter(class != "all_genes") %>%
  ggplot(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),prop.sel,fill=factor(class,levels=c("all_genes","monophyletic","deep_paraphyletic","shallow_paraphyletic")))) +
  geom_bar(stat = "identity",position=position_dodge()) +
  geom_text(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),label=sig),position=position_dodge(width=0.9),size=10) +
  scale_x_discrete(labels=c("0.1","0.01","0.001","0.0001")) +
  scale_fill_discrete(name="class") +
  xlab("q-value") +
  ylab("proportion selected birds + mammals") +
  labs(subtitle="* = p < 0.05, ** = p < 0.01, *** = p < 0.001") +
  facet_grid(~factor(comparison,levels=c("mammals","birds","both")))

odds_ratio_plot <- qval_res_trans_forplotting %>%
  filter(class != "all_genes") %>%
  ggplot(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),odds.ratio,col=factor(class,levels=c("monophyletic","deep_paraphyletic","shallow_paraphyletic")))) +
  geom_point(size=4,position=position_dodge(width=0.9)) +
  geom_linerange(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),ymin=conf.lower,ymax=conf.upper,col=factor(class,levels=c("monophyletic","deep_paraphyletic","shallow_paraphyletic"))),position=position_dodge(width=0.9),size=1.5) +
  geom_hline(aes(yintercept = 1)) +
  scale_x_discrete(labels=c("0.1","0.01","0.001","0.0001")) +
  scale_color_discrete(name="class") +
  xlab("q-value") +
  ylab("odds ratio") +
  scale_y_continuous(labels=function(label) sprintf('%6.0f', label)) +
  facet_grid(~factor(comparison,levels=c("mammals","birds","both")))

plot_grid(prop_sel_plot, odds_ratio_plot, ncol=1)
ggsave(filename = "07_output_virus_classification//mammal_bird_evidence_type_comparison.pdf",width = 14, height=10)


#######################################################################################################################
##################### Second, look for enrichment in different viral families #####################
#######################################################################################################################

#q values to consider
qvals <- c(0.1,0.01,0.001,0.0001)

#List to capture results across q values
qval_res_family_list <- list()

for (i in 1:length(qvals)){
  
  comp_propsig_family <- matrix(nrow=(length(families_to_test)+1),ncol=18)
  comp_propsig_family[,1] <- c(families_to_test,"all_genes")
  
  for (j in 1:length(families_to_test)){
    #Create a vector of true and false for that evidence type for each hog, and remove duplicates
    family_hogs = imm %>%
      filter(family == families_to_test[j]) %>%
      pull(hog)
    imm_test <- imm %>%
      distinct(hog,.keep_all=TRUE) %>%
      mutate(family_type = if_else(hog %in% family_hogs,TRUE,FALSE))
    
    #Calculate fisher's exact test for prop selected in both birds and mammals, and get prop sig genes
    comp_propsig_family[j,2:5] <- imm_test %>%
      with(.,table(family_type,mammal_q < qvals[i] & bird_q < qvals[i])) %>% fisher.test %>% unlist %>% .[1:4]
    comp_propsig_family[j,6] <- imm_test %>%
      filter(family_type == TRUE) %>%
      with(.,prop.table(table(mammal_q < qvals[i] & bird_q < qvals[i]))) %>% .[2]
    
    #Calculate fisher's exact test for prop selected in mammals only, and get prop sig genes
    comp_propsig_family[j,7:10] <- imm_test %>%
      with(.,table(family_type,mammal_q < qvals[i])) %>% fisher.test %>% unlist %>% .[1:4]
    comp_propsig_family[j,11] <- imm_test %>%
      filter(family_type == TRUE) %>%
      with(.,prop.table(table(mammal_q < qvals[i]))) %>% .[2]
    
    #Calculate fisher's exact test for prop selected in birds only, and get prop sig genes
    comp_propsig_family[j,12:15] <- imm_test %>%
      with(.,table(family_type,bird_q < qvals[i])) %>% fisher.test %>% unlist %>% .[1:4]
    comp_propsig_family[j,16] <- imm_test %>%
      filter(family_type == TRUE) %>%
      with(.,prop.table(table(bird_q < qvals[i]))) %>% .[2]
    
    #Number of genes in each category:
    comp_propsig_family[j,17] <- imm_test %>% filter(family_type == TRUE) %>% summarize(n()) %>% pull
  }
  
  #Add prop sig all genes
  comp_propsig_family[j+1,6] <- imm_test %>%
    with(.,prop.table(table(mammal_q < qvals[i] & bird_q < qvals[i]))) %>% .[2]
  comp_propsig_family[j+1,11] <- imm_test %>%
    with(.,prop.table(table(mammal_q < qvals[i]))) %>% .[2]
  comp_propsig_family[j+1,16] <- imm_test %>%
    with(.,prop.table(table(bird_q < qvals[i]))) %>% .[2]
  comp_propsig_family[j+1,17] <- imm_test %>% summarize(n()) %>% pull
  
  comp_propsig_family[,18] <- qvals[i]
  
  colnames(comp_propsig_family) <- c("class","p.value_both","conf.int1_both","conf.int2_both","estimated.odds.ratio_both","prop_sel_both","p.value_mammals","conf.int1_mammals","conf.int2_mammals","estimated.odds.ratio_mammals","prop_sel_mammals","p.value_birds","conf.int1_birds","conf.int2_birds","estimated.odds.ratio_birds","prop_sel_birds","n.genes","qval")
  
  #Clean up, select relevant columns
  comp_propsig_family_clean <- comp_propsig_family %>%
    as.tibble %>%
    mutate(p.value_both=round(as.numeric(p.value_both),digits = 4),
           odds.ratio_both=round(as.double(estimated.odds.ratio_both),digits=2),
           prop.sel_both=round(as.numeric(prop_sel_both),digits=3),
           conf.lower_both=round(as.numeric(conf.int1_both),digits=3),
           conf.upper_both=round(as.numeric(conf.int2_both),digits=3),
           p.value_mammals=round(as.numeric(p.value_mammals),digits = 4),
           odds.ratio_mammals=round(as.numeric(estimated.odds.ratio_mammals),digits=2),
           prop.sel_mammals=round(as.numeric(prop_sel_mammals),digits=3),
           conf.lower_mammals=round(as.numeric(conf.int1_mammals),digits=3),
           conf.upper_mammals=round(as.numeric(conf.int2_mammals),digits=3),
           p.value_birds=round(as.numeric(p.value_birds),digits = 4),
           odds.ratio_birds=round(as.numeric(estimated.odds.ratio_birds),digits=2),
           prop.sel_birds=round(as.numeric(prop_sel_birds),digits=3),
           conf.lower_birds=round(as.numeric(conf.int1_birds),digits=3),
           conf.upper_birds=round(as.numeric(conf.int2_birds),digits=3)) %>%
    dplyr::select(qval,class,n.genes,odds.ratio_both,conf.lower_both,conf.upper_both,p.value_both,prop.sel_both,odds.ratio_mammals,conf.lower_mammals,conf.upper_mammals,p.value_mammals,prop.sel_mammals,odds.ratio_birds,conf.lower_birds,conf.upper_birds,p.value_birds,prop.sel_birds)
  
  
  qval_res_family_list[[i]] <- comp_propsig_family_clean
}

qval_res_family <- qval_res_family_list %>% bind_rows

write_csv(qval_res_family,path="07_output_virus_classification/bird_mammals_viral_family_enrichment_table.csv")


#Reshape into long format for plotting
qval_res_family <- qval_res_family %>%
  gather(stat,value,-qval,-class,-n.genes) %>%
  separate(stat,into = c("stat","comparison"),sep="_")

#Create vector of asterisks to include in plots:
qval_res_family_forplotting <-   qval_res_family %>%
  #filter(class != "all_genes") %>%
  spread(stat,value) %>%
  mutate(sig = case_when(
    is.na(p.value) ~ "",
    p.value > 0.05 ~ "",
    p.value <= 0.05 & p.value > 0.01 ~ "*",
    p.value <= 0.01 & p.value > 0.001 ~ "**",
    p.value <= 0.001 ~ "***"
  ))


prop_sel_plot <- qval_res_family_forplotting %>%
  filter(class != "all_genes") %>%
  ggplot(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),prop.sel,fill=factor(class,levels=c("all_genes",families_to_test)))) +
  geom_bar(stat = "identity",position=position_dodge()) +
  geom_text(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),label=sig),position=position_dodge(width=0.9),size=10) +
  scale_x_discrete(labels=c("0.1","0.01","0.001","0.0001")) +
  scale_fill_discrete(name="class") +
  xlab("q-value") +
  ylab("proportion selected birds + mammals") +
  labs(subtitle="* = p < 0.05, ** = p < 0.01, *** = p < 0.001") +
  facet_grid(~factor(comparison,levels=c("mammals","birds","both")))

odds_ratio_plot <- qval_res_family_forplotting %>%
  filter(class != "all_genes") %>%
  ggplot(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),odds.ratio,col=factor(class,levels=families_to_test))) +
  geom_point(size=4,position=position_dodge(width=0.9)) +
  geom_linerange(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),ymin=conf.lower,ymax=conf.upper,col=factor(class,levels=families_to_test)),position=position_dodge(width=0.9),size=1.5) +
  geom_hline(aes(yintercept = 1)) +
  scale_x_discrete(labels=c("0.1","0.01","0.001","0.0001")) +
  scale_color_discrete(name="class") +
  xlab("q-value") +
  ylab("odds ratio") +
  scale_y_continuous(labels=function(label) sprintf('%6.0f', label)) +
  facet_grid(~factor(comparison,levels=c("mammals","birds","both")))

plot_grid(prop_sel_plot, odds_ratio_plot, ncol=1)
ggsave(filename = "07_output_virus_classification/mammal_bird_family_comparison.pdf",width = 18, height=10)

#Plot with joint results only
prop_sel_plot <- qval_res_family_forplotting %>%
  filter(class != "all_genes") %>%
  filter(comparison == "both") %>%
  ggplot(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),prop.sel,fill=factor(class,levels=c("all_genes",families_to_test)))) +
  geom_bar(stat = "identity",position=position_dodge()) +
  geom_text(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),label=sig),position=position_dodge(width=0.9),size=10) +
  scale_x_discrete(labels=c("0.1","0.01","0.001","0.0001")) +
  scale_fill_discrete(name="class") +
  xlab("q-value") +
  ylab("proportion selected birds + mammals") +
  labs(subtitle="* = p < 0.05, ** = p < 0.01, *** = p < 0.001")

odds_ratio_plot <- qval_res_family_forplotting %>%
  filter(class != "all_genes") %>%
  filter(comparison == "both") %>%
  ggplot(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),odds.ratio,col=factor(class,levels=families_to_test))) +
  geom_point(size=4,position=position_dodge(width=0.9)) +
  geom_linerange(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),ymin=conf.lower,ymax=conf.upper,col=factor(class,levels=families_to_test)),position=position_dodge(width=0.9),size=1.5) +
  geom_hline(aes(yintercept = 1)) +
  scale_x_discrete(labels=c("0.1","0.01","0.001","0.0001")) +
  scale_color_discrete(name="class") +
  xlab("q-value") +
  ylab("odds ratio") +
  scale_y_continuous(labels=function(label) sprintf('%6.0f', label))

plot_grid(prop_sel_plot, odds_ratio_plot, ncol=1)
ggsave(filename = "07_output_virus_classification/mammal_bird_family_comparison_joint_results_only.pdf",width = 10, height=10)



#######################################################################################################################
##################### Third, look for enrichment in different classes of evidence types compared to VIPs #####################
#######################################################################################################################


#q values to consider
qvals <- c(0.1,0.01,0.001,0.0001)

#List to capture results across q values
qval_res_trans_list <- list()

for (i in 1:length(qvals)){
  
  comp_propsig_trans <- matrix(nrow=(length(trans_evidence_types)+1),ncol=18)
  comp_propsig_trans[,1] <- c(trans_evidence_types,"VIPs")
  
  for (j in 1:length(trans_evidence_types)){
    #Create a vector of true and false for that evidence type for each hog, and remove duplicates
    trans_evi_type_hogs = imm %>%
      filter(trans_evidence == trans_evidence_types[j]) %>%
      pull(hog)
    imm_test <- imm %>%
      distinct(hog,.keep_all=TRUE) %>%
      filter(vip) %>%
      mutate(trans_type = if_else(hog %in% trans_evi_type_hogs,TRUE,FALSE))
    
    #Calculate fisher's exact test for prop selected in both birds and mammals, and get prop sig genes
    comp_propsig_trans[j,2:5] <- imm_test %>%
      with(.,table(trans_type,mammal_q < qvals[i] & bird_q < qvals[i])) %>% fisher.test %>% unlist %>% .[1:4]
    comp_propsig_trans[j,6] <- imm_test %>%
      filter(trans_type == TRUE) %>%
      with(.,prop.table(table(mammal_q < qvals[i] & bird_q < qvals[i]))) %>% .[2]
    
    #Calculate fisher's exact test for prop selected in mammals only, and get prop sig genes
    comp_propsig_trans[j,7:10] <- imm_test %>%
      with(.,table(trans_type,mammal_q < qvals[i])) %>% fisher.test %>% unlist %>% .[1:4]
    comp_propsig_trans[j,11] <- imm_test %>%
      filter(trans_type == TRUE) %>%
      with(.,prop.table(table(mammal_q < qvals[i]))) %>% .[2]
    
    #Calculate fisher's exact test for prop selected in birds only, and get prop sig genes
    comp_propsig_trans[j,12:15] <- imm_test %>%
      with(.,table(trans_type,bird_q < qvals[i])) %>% fisher.test %>% unlist %>% .[1:4]
    comp_propsig_trans[j,16] <- imm_test %>%
      filter(trans_type == TRUE) %>%
      with(.,prop.table(table(bird_q < qvals[i]))) %>% .[2]
    
    #Number of genes in each category:
    comp_propsig_trans[j,17] <- imm_test %>% filter(trans_type == TRUE) %>% summarize(n()) %>% pull
  }
  
  #Add prop sig all genes
  comp_propsig_trans[j+1,6] <- imm_test %>%
    with(.,prop.table(table(mammal_q < qvals[i] & bird_q < qvals[i]))) %>% .[2]
  comp_propsig_trans[j+1,11] <- imm_test %>%
    with(.,prop.table(table(mammal_q < qvals[i]))) %>% .[2]
  comp_propsig_trans[j+1,16] <- imm_test %>%
    with(.,prop.table(table(bird_q < qvals[i]))) %>% .[2]
  comp_propsig_trans[j+1,17] <- imm_test %>% summarize(n()) %>% pull
  
  comp_propsig_trans[,18] <- qvals[i]
  
  colnames(comp_propsig_trans) <- c("class","p.value_both","conf.int1_both","conf.int2_both","estimated.odds.ratio_both","prop_sel_both","p.value_mammals","conf.int1_mammals","conf.int2_mammals","estimated.odds.ratio_mammals","prop_sel_mammals","p.value_birds","conf.int1_birds","conf.int2_birds","estimated.odds.ratio_birds","prop_sel_birds","n.genes","qval")
  
  #Clean up, select relevant columns
  comp_propsig_trans_clean <- comp_propsig_trans %>%
    as.tibble %>%
    mutate(p.value_both=round(as.numeric(p.value_both),digits = 4),
           odds.ratio_both=round(as.double(estimated.odds.ratio_both),digits=2),
           prop.sel_both=round(as.numeric(prop_sel_both),digits=3),
           conf.lower_both=round(as.numeric(conf.int1_both),digits=3),
           conf.upper_both=round(as.numeric(conf.int2_both),digits=3),
           p.value_mammals=round(as.numeric(p.value_mammals),digits = 4),
           odds.ratio_mammals=round(as.numeric(estimated.odds.ratio_mammals),digits=2),
           prop.sel_mammals=round(as.numeric(prop_sel_mammals),digits=3),
           conf.lower_mammals=round(as.numeric(conf.int1_mammals),digits=3),
           conf.upper_mammals=round(as.numeric(conf.int2_mammals),digits=3),
           p.value_birds=round(as.numeric(p.value_birds),digits = 4),
           odds.ratio_birds=round(as.numeric(estimated.odds.ratio_birds),digits=2),
           prop.sel_birds=round(as.numeric(prop_sel_birds),digits=3),
           conf.lower_birds=round(as.numeric(conf.int1_birds),digits=3),
           conf.upper_birds=round(as.numeric(conf.int2_birds),digits=3)) %>%
    dplyr::select(qval,class,n.genes,odds.ratio_both,conf.lower_both,conf.upper_both,p.value_both,prop.sel_both,odds.ratio_mammals,conf.lower_mammals,conf.upper_mammals,p.value_mammals,prop.sel_mammals,odds.ratio_birds,conf.lower_birds,conf.upper_birds,p.value_birds,prop.sel_birds)
  
  
  qval_res_trans_list[[i]] <- comp_propsig_trans_clean
}

qval_res_trans <- qval_res_trans_list %>% bind_rows

write_csv(qval_res_trans,path="07_output_virus_classification/bird_mammals_evidence_type_enrichment_table_vips_only.csv")


#Reshape into long format for plotting
qval_res_trans <- qval_res_trans %>%
  gather(stat,value,-qval,-class,-n.genes) %>%
  separate(stat,into = c("stat","comparison"),sep="_")

#Create vector of asterisks to include in plots:
qval_res_trans_forplotting <-   qval_res_trans %>%
  #filter(class != "all_genes") %>%
  filter(class != "only_mammals") %>%
  spread(stat,value) %>%
  mutate(sig = case_when(
    is.na(p.value) ~ "",
    p.value > 0.05 ~ "",
    p.value <= 0.05 & p.value > 0.01 ~ "*",
    p.value <= 0.01 & p.value > 0.001 ~ "**",
    p.value <= 0.001 ~ "***"
  ))


prop_sel_plot <- qval_res_trans_forplotting %>%
  filter(class != "VIPs") %>%
  ggplot(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),prop.sel,fill=factor(class,levels=c("all_genes","monophyletic","deep_paraphyletic","shallow_paraphyletic")))) +
  geom_bar(stat = "identity",position=position_dodge()) +
  geom_text(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),label=sig),position=position_dodge(width=0.9),size=10) +
  scale_x_discrete(labels=c("0.1","0.01","0.001","0.0001")) +
  scale_fill_discrete(name="class") +
  xlab("q-value") +
  ylab("proportion selected birds + mammals") +
  labs(subtitle="* = p < 0.05, ** = p < 0.01, *** = p < 0.001") +
  facet_grid(~factor(comparison,levels=c("mammals","birds","both")))

odds_ratio_plot <- qval_res_trans_forplotting %>%
  filter(class != "VIPs") %>%
  ggplot(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),odds.ratio,col=factor(class,levels=c("monophyletic","deep_paraphyletic","shallow_paraphyletic")))) +
  geom_point(size=4,position=position_dodge(width=0.9)) +
  geom_linerange(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),ymin=conf.lower,ymax=conf.upper,col=factor(class,levels=c("monophyletic","deep_paraphyletic","shallow_paraphyletic"))),position=position_dodge(width=0.9),size=1.5) +
  geom_hline(aes(yintercept = 1)) +
  scale_x_discrete(labels=c("0.1","0.01","0.001","0.0001")) +
  scale_color_discrete(name="class") +
  xlab("q-value") +
  ylab("odds ratio") +
  scale_y_continuous(labels=function(label) sprintf('%6.0f', label)) +
  facet_grid(~factor(comparison,levels=c("mammals","birds","both")))

plot_grid(prop_sel_plot, odds_ratio_plot, ncol=1)
ggsave(filename = "07_output_virus_classification//mammal_bird_evidence_type_comparison_vips_only.pdf",width = 14, height=10)


#######################################################################################################################
##################### Second, look for enrichment in different viral families compared to VIPs #####################
#######################################################################################################################

#q values to consider
qvals <- c(0.1,0.01,0.001,0.0001)

#List to capture results across q values
qval_res_family_list <- list()

for (i in 1:length(qvals)){
  
  comp_propsig_family <- matrix(nrow=(length(families_to_test)+1),ncol=18)
  comp_propsig_family[,1] <- c(families_to_test,"VIPs")
  
  for (j in 1:length(families_to_test)){
    #Create a vector of true and false for that evidence type for each hog, and remove duplicates
    family_hogs = imm %>%
      filter(family == families_to_test[j]) %>%
      pull(hog)
    imm_test <- imm %>%
      distinct(hog,.keep_all=TRUE) %>%
      filter(vip) %>%
      mutate(family_type = if_else(hog %in% family_hogs,TRUE,FALSE))
    
    #Calculate fisher's exact test for prop selected in both birds and mammals, and get prop sig genes
    comp_propsig_family[j,2:5] <- imm_test %>%
      with(.,table(family_type,mammal_q < qvals[i] & bird_q < qvals[i])) %>% fisher.test %>% unlist %>% .[1:4]
    comp_propsig_family[j,6] <- imm_test %>%
      filter(family_type == TRUE) %>%
      with(.,prop.table(table(mammal_q < qvals[i] & bird_q < qvals[i]))) %>% .[2]
    
    #Calculate fisher's exact test for prop selected in mammals only, and get prop sig genes
    comp_propsig_family[j,7:10] <- imm_test %>%
      with(.,table(family_type,mammal_q < qvals[i])) %>% fisher.test %>% unlist %>% .[1:4]
    comp_propsig_family[j,11] <- imm_test %>%
      filter(family_type == TRUE) %>%
      with(.,prop.table(table(mammal_q < qvals[i]))) %>% .[2]
    
    #Calculate fisher's exact test for prop selected in birds only, and get prop sig genes
    comp_propsig_family[j,12:15] <- imm_test %>%
      with(.,table(family_type,bird_q < qvals[i])) %>% fisher.test %>% unlist %>% .[1:4]
    comp_propsig_family[j,16] <- imm_test %>%
      filter(family_type == TRUE) %>%
      with(.,prop.table(table(bird_q < qvals[i]))) %>% .[2]
    
    #Number of genes in each category:
    comp_propsig_family[j,17] <- imm_test %>% filter(family_type == TRUE) %>% summarize(n()) %>% pull
  }
  
  #Add prop sig all genes
  comp_propsig_family[j+1,6] <- imm_test %>%
    with(.,prop.table(table(mammal_q < qvals[i] & bird_q < qvals[i]))) %>% .[2]
  comp_propsig_family[j+1,11] <- imm_test %>%
    with(.,prop.table(table(mammal_q < qvals[i]))) %>% .[2]
  comp_propsig_family[j+1,16] <- imm_test %>%
    with(.,prop.table(table(bird_q < qvals[i]))) %>% .[2]
  comp_propsig_family[j+1,17] <- imm_test %>% summarize(n()) %>% pull
  
  comp_propsig_family[,18] <- qvals[i]
  
  colnames(comp_propsig_family) <- c("class","p.value_both","conf.int1_both","conf.int2_both","estimated.odds.ratio_both","prop_sel_both","p.value_mammals","conf.int1_mammals","conf.int2_mammals","estimated.odds.ratio_mammals","prop_sel_mammals","p.value_birds","conf.int1_birds","conf.int2_birds","estimated.odds.ratio_birds","prop_sel_birds","n.genes","qval")
  
  #Clean up, select relevant columns
  comp_propsig_family_clean <- comp_propsig_family %>%
    as.tibble %>%
    mutate(p.value_both=round(as.numeric(p.value_both),digits = 4),
           odds.ratio_both=round(as.double(estimated.odds.ratio_both),digits=2),
           prop.sel_both=round(as.numeric(prop_sel_both),digits=3),
           conf.lower_both=round(as.numeric(conf.int1_both),digits=3),
           conf.upper_both=round(as.numeric(conf.int2_both),digits=3),
           p.value_mammals=round(as.numeric(p.value_mammals),digits = 4),
           odds.ratio_mammals=round(as.numeric(estimated.odds.ratio_mammals),digits=2),
           prop.sel_mammals=round(as.numeric(prop_sel_mammals),digits=3),
           conf.lower_mammals=round(as.numeric(conf.int1_mammals),digits=3),
           conf.upper_mammals=round(as.numeric(conf.int2_mammals),digits=3),
           p.value_birds=round(as.numeric(p.value_birds),digits = 4),
           odds.ratio_birds=round(as.numeric(estimated.odds.ratio_birds),digits=2),
           prop.sel_birds=round(as.numeric(prop_sel_birds),digits=3),
           conf.lower_birds=round(as.numeric(conf.int1_birds),digits=3),
           conf.upper_birds=round(as.numeric(conf.int2_birds),digits=3)) %>%
    dplyr::select(qval,class,n.genes,odds.ratio_both,conf.lower_both,conf.upper_both,p.value_both,prop.sel_both,odds.ratio_mammals,conf.lower_mammals,conf.upper_mammals,p.value_mammals,prop.sel_mammals,odds.ratio_birds,conf.lower_birds,conf.upper_birds,p.value_birds,prop.sel_birds)
  
  
  qval_res_family_list[[i]] <- comp_propsig_family_clean
}

qval_res_family <- qval_res_family_list %>% bind_rows

write_csv(qval_res_family,path="07_output_virus_classification/bird_mammals_viral_family_enrichment_table_vips_only.csv")


#Reshape into long format for plotting
qval_res_family <- qval_res_family %>%
  gather(stat,value,-qval,-class,-n.genes) %>%
  separate(stat,into = c("stat","comparison"),sep="_")

#Create vector of asterisks to include in plots:
qval_res_family_forplotting <-   qval_res_family %>%
  #filter(class != "all_genes") %>%
  spread(stat,value) %>%
  mutate(sig = case_when(
    is.na(p.value) ~ "",
    p.value > 0.05 ~ "",
    p.value <= 0.05 & p.value > 0.01 ~ "*",
    p.value <= 0.01 & p.value > 0.001 ~ "**",
    p.value <= 0.001 ~ "***"
  ))


prop_sel_plot <- qval_res_family_forplotting %>%
  #filter(class != "VIPs") %>%
  ggplot(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),prop.sel,fill=factor(class,levels=c("VIPs",families_to_test)))) +
  geom_bar(stat = "identity",position=position_dodge()) +
  geom_text(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),label=sig),position=position_dodge(width=0.9),size=10) +
  scale_x_discrete(labels=c("0.1","0.01","0.001","0.0001")) +
  scale_fill_discrete(name="class") +
  xlab("q-value") +
  ylab("proportion selected birds + mammals") +
  labs(subtitle="* = p < 0.05, ** = p < 0.01, *** = p < 0.001") +
  facet_grid(~factor(comparison,levels=c("mammals","birds","both")))

odds_ratio_plot <- qval_res_family_forplotting %>%
  filter(class != "VIPs") %>%
  ggplot(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),odds.ratio,col=factor(class,levels=families_to_test))) +
  geom_point(size=4,position=position_dodge(width=0.9)) +
  geom_linerange(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),ymin=conf.lower,ymax=conf.upper,col=factor(class,levels=families_to_test)),position=position_dodge(width=0.9),size=1.5) +
  geom_hline(aes(yintercept = 1)) +
  scale_x_discrete(labels=c("0.1","0.01","0.001","0.0001")) +
  scale_color_discrete(name="class") +
  xlab("q-value") +
  ylab("odds ratio") +
  scale_y_continuous(labels=function(label) sprintf('%6.0f', label)) +
  facet_grid(~factor(comparison,levels=c("mammals","birds","both")))

plot_grid(prop_sel_plot, odds_ratio_plot, ncol=1)
ggsave(filename = "07_output_virus_classification/mammal_bird_family_comparison_vips_only.pdf",width = 18, height=10)

#Plot with joint results only
prop_sel_plot <- qval_res_family_forplotting %>%
  #filter(class != "VIPs") %>%
  filter(comparison == "both") %>%
  ggplot(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),prop.sel,fill=factor(class,levels=c("VIPs",families_to_test)))) +
  geom_bar(stat = "identity",position=position_dodge()) +
  geom_text(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),label=sig),position=position_dodge(width=0.9),size=10) +
  scale_x_discrete(labels=c("0.1","0.01","0.001","0.0001")) +
  scale_fill_discrete(name="class") +
  xlab("q-value") +
  ylab("proportion selected birds + mammals") +
  labs(subtitle="* = p < 0.05, ** = p < 0.01, *** = p < 0.001")

odds_ratio_plot <- qval_res_family_forplotting %>%
  filter(class != "VIPs") %>%
  filter(comparison == "both") %>%
  ggplot(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),odds.ratio,col=factor(class,levels=families_to_test))) +
  geom_point(size=4,position=position_dodge(width=0.9)) +
  geom_linerange(aes(factor(qval,levels=c("0.1","0.01","0.001","1e-04")),ymin=conf.lower,ymax=conf.upper,col=factor(class,levels=families_to_test)),position=position_dodge(width=0.9),size=1.5) +
  geom_hline(aes(yintercept = 1)) +
  scale_x_discrete(labels=c("0.1","0.01","0.001","0.0001")) +
  scale_color_discrete(name="class") +
  xlab("q-value") +
  ylab("odds ratio") +
  scale_y_continuous(labels=function(label) sprintf('%6.0f', label))

plot_grid(prop_sel_plot, odds_ratio_plot, ncol=1)
ggsave(filename = "07_output_virus_classification/mammal_bird_family_comparison_joint_results_only_vips_only.pdf",width = 10, height=10)

