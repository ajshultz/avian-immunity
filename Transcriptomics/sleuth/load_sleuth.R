#Read in sleuth data
#Loosely based on https://github.com/pachterlab/bears_analyses/blob/master/sleuth.R

#Preparation (load libraries and target id to gene id data)
library(tidyverse)
library(sleuth)
t2g<-read_tsv("../t2g_key.clean", col_names = c("target_id", "gene_id"))

#Function to prep data based on a BioProject ID

load_samples_for_sleuth <- function(bioproject) {
  print(paste0("Working on BioProject: ", bioproject, collapse=""))
  file <- paste0(bioproject, ".prep", collapse="")
  read_tsv(file) %>% distinct
}

#Function to fit full and reduced model to data

load_sleuth_obj <- function(s2c) {
  so<-sleuth_prep(s2c, target_mapping=t2g, num_cores=4, aggregation_column = "gene_id", extra_bootstrap_summary = TRUE, transformation_function = function(x) log2(x + 0.5))
  return(so)
}

export_pca <- function(so, condition, bioproject) {
  plot_pca(so, color_by = condition, units = "scaled_reads_per_base", text_labels = TRUE)
  file <- paste0(bioproject, "-", condition, ".pca.pdf", collapse="")
  sampfile <- paste0(bioproject, ".s2c.tsv", collapse="")
  ggsave(file)
  so$sample_to_covariates %>% write_tsv(sampfile)
}

export_results <- function(so, condition, model, bioproject, infect) {
  file <- paste0(bioproject, "-", condition, ".results.tsv", collapse="")
  vplot <- paste0(bioproject, "-", condition, ".volcano.pdf", collapse="")
  so<-sleuth_wt(so, condition, model)
  sleuth_results(so, condition, which_model = model, show_all=TRUE) %>% 
    select(target_id:var_obs) %>%
    mutate(bioproject = bioproject, condition = condition, infect = infect) %>% 
    mutate(sig = sign(b) * (qval <= 0.05)) %>%
    mutate(species = ifelse(substr(target_id, 1,3) == "ENS", 
                            ifelse(substr(target_id, 4, 6) == "G00", "HOM", substr(target_id, 4,6)), "SCA")) %>%
    write_tsv(file)
  plot_volcano(so, condition, which_model = model, sig_level = 0.05)
  ggsave(vplot)
}

#Steps for each BioProject
#1: load all samples
#2: split into subsets for testing
#3: run sleuth on subsets
#4: export results and plots

##BIOPROJECT: PRJEB1284

bioproject<-"PRJEB1284"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$timepoint = as.factor(allsamples$timepoint)
allsamples$timepoint = relevel(allsamples$timepoint, ref="0")
allsamples$treatment = as.factor(allsamples$treatment)
allsamples$treatment = relevel(allsamples$treatment, ref="Control")
so<-load_sleuth_obj(allsamples)

#pca
export_pca(so, "treatment", bioproject)
export_pca(so, "timepoint", bioproject)

#analysis
#treatment plas vs control
so<-sleuth_fit(so, ~treatment, "treatment")
export_results(so, "treatmentPlasmodium_chabaudi", "treatment", bioproject, "Plasmodium_chabaudi")

#timepoints vs 0
so<-sleuth_fit(so, ~timepoint, "timepoint")
export_results(so, "timepoint6", "timepoint", bioproject, "Plasmodium_chabaudi")
export_results(so, "timepoint7", "timepoint", bioproject, "Plasmodium_chabaudi")
export_results(so, "timepoint8", "timepoint", bioproject, "Plasmodium_chabaudi")
export_results(so, "timepoint9", "timepoint", bioproject, "Plasmodium_chabaudi")
export_results(so, "timepoint10", "timepoint", bioproject, "Plasmodium_chabaudi")
export_results(so, "timepoint11", "timepoint", bioproject, "Plasmodium_chabaudi")
export_results(so, "timepoint12", "timepoint", bioproject, "Plasmodium_chabaudi")

###################

##BIOPROJECT: PRJEB4731

bioproject<-"PRJEB4731"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$timepoint = as.factor(allsamples$timepoint)
allsamples$timepoint = relevel(allsamples$timepoint, ref="0")

#fix issues -- C3H and C3H/He are same, control capitalization not consistent
allsamples <- allsamples %>% mutate(treatment = ifelse(treatment=="control", "Control", treatment)) %>%
  mutate(genotype = ifelse(genotype=="C3H/He", "C3H", genotype)) %>%
  mutate(treatment = relevel(as.factor(treatment), ref="Control"))


#will use genotype as a block variable by default but check PCA
so<-load_sleuth_obj(allsamples)

#pca
export_pca(so, "treatment", bioproject)
export_pca(so, "timepoint", bioproject)
export_pca(so, "genotype", bioproject)

#analysis
#treatment plas vs control
so<-sleuth_fit(so, ~treatment+genotype, "treatment")
export_results(so, "treatmentPlasmodium_chabaudi", "treatment", bioproject, "Plasmodium_chabaudi")

#timepoints vs 0
so<-sleuth_fit(so, ~timepoint+genotype, "timepoint")
export_results(so, "timepoint3", "timepoint", bioproject, "Plasmodium_chabaudi")
export_results(so, "timepoint5", "timepoint", bioproject, "Plasmodium_chabaudi")
export_results(so, "timepoint6", "timepoint", bioproject, "Plasmodium_chabaudi")
export_results(so, "timepoint7", "timepoint", bioproject, "Plasmodium_chabaudi")
export_results(so, "timepoint8", "timepoint", bioproject, "Plasmodium_chabaudi")
export_results(so, "timepoint9", "timepoint", bioproject, "Plasmodium_chabaudi")
export_results(so, "timepoint10", "timepoint", bioproject, "Plasmodium_chabaudi")
export_results(so, "timepoint11", "timepoint", bioproject, "Plasmodium_chabaudi")
export_results(so, "timepoint12", "timepoint", bioproject, "Plasmodium_chabaudi")


###################

##BIOPROJECT: PRJEB7213

bioproject<-"PRJEB7213"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$treatment = relevel(as.factor(allsamples$treatment), ref="none")

#split into day/tissue, so four experiments

for (i in c(1, 3)) {
  for (j in c("lung", "ileum")) {
    bioproject2 <- paste0(bioproject, "-time", i, "-", j, collapse = "")
    subsample <- allsamples %>% filter(time == i, tissue == j)
    so<-load_sleuth_obj(subsample)
    export_pca(so, "treatment", bioproject2)
    so<-sleuth_fit(so, ~treatment, "treatment")
    export_results(so, "treatmentH5N2", "treatment", bioproject2, "H5N2")
    export_results(so, "treatmentH5N1", "treatment", bioproject2, "H5N1")
  }
}

###################

##BIOPROJECT: PRJEB7215


bioproject<-"PRJEB7215"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$treatment = relevel(as.factor(allsamples$treatment), ref="none")

#split into day/tissue, so four experiments

for (i in c(1, 3)) {
  for (j in c("lung", "ileum")) {
    bioproject2 <- paste0(bioproject, "-time", i, "-", j, collapse = "")
    subsample <- allsamples %>% filter(time == i, tissue == j)
    so<-load_sleuth_obj(subsample)
    export_pca(so, "treatment", bioproject2)
    so<-sleuth_fit(so, ~treatment, "treatment")
    export_results(so, "treatmentH5N2", "treatment", bioproject2, "H5N2")
    export_results(so, "treatmentH5N1", "treatment", bioproject2, "H5N1")
  }
}

###################

##BIOPROJECT: PRJEB7219

bioproject<-"PRJEB7219"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$treatment = relevel(as.factor(allsamples$treatment), ref="none")

so<-load_sleuth_obj(allsamples)
export_pca(so, "treatment", bioproject)
so<-sleuth_fit(so, ~treatment, "treatment")
export_results(so, "treatmentIBDV", "treatment", bioproject, "IBDV")

###################

##BIOPROJECT: PRJNA227074

bioproject<-"PRJNA227074"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$treatment = relevel(as.factor(allsamples$treatment), ref="none")

#three outliers in PCA, looks driven by mtDNA mappings
#SAMN02665938 - ind 210
#SAMN02665937 - ind 212
#SAMN02665945 - ind 233

#remove these individuals
all_sub <- allsamples %>% filter(individual == 210 | individual == 212 | individual == 233)
so<-load_sleuth_obj(all_sub)
export_pca(so, "treatment", bioproject)
so<-sleuth_fit(so, ~treatment, "treatment")
export_results(so, "treatmentPlasmodium_falciparum", "treatment", bioproject, "Plasmodium_falciparum")

###################

##BIOPROJECT: PRJEB17676

bioproject<-"PRJEB17676"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$treatment = relevel(as.factor(allsamples$treatment), ref="Control")

#split by tissue
for (i in c("mammary_gland", "spleen")) {
  bioproject2 <- paste0(bioproject, "-", i, collapse = "")
  subsample <- allsamples %>% filter(tissue == i)
  so<-load_sleuth_obj(subsample)
  export_pca(so, "treatment", bioproject2)
  so<-sleuth_fit(so, ~treatment, "treatment")
  export_results(so, "treatmentMycoplasma_agalactiae", "treatment", bioproject2, "Mycoplasma_agalactiae")
}

###################


##BIOPROJECT: PRJEB19318

bioproject<-"PRJEB19318"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$treatment = relevel(as.factor(allsamples$treatment), ref="control")

#split by genotype and time, treat sex as blocking variable
for (i in c(2, 6, 10)) {
  for (j in c("resistant", "susceptible")) {
    bioproject2 <- paste0(bioproject, "-time", i, "-", j, collapse = "")
    subsample <- allsamples %>% filter(time == i, genotype == j)
    so<-load_sleuth_obj(subsample)
    export_pca(so, "treatment", bioproject2)
    export_pca(so, "sex", bioproject2)
    so<-sleuth_fit(so, ~treatment+sex, "treatment")
    export_results(so, "treatmentNDV", "treatment", bioproject2, "NDV")
  }
}

###################

##BIOPROJECT: PRJEB21688

bioproject<-"PRJEB21688"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$treatment = relevel(as.factor(allsamples$treatment), ref="control")

#split by genotype and time, treat sex as blocking variable
for (i in c(2, 6)) {
  for (j in c("resistant", "susceptible")) {
    bioproject2 <- paste0(bioproject, "-time", i, "-", j, collapse = "")
    subsample <- allsamples %>% filter(time == i, genotype == j)
    so<-load_sleuth_obj(subsample)
    export_pca(so, "treatment", bioproject2)
    export_pca(so, "sex", bioproject2)
    so<-sleuth_fit(so, ~treatment+sex, "treatment")
    export_results(so, "treatmentNDV", "treatment", bioproject2, "NDV")
  }
}

###################

##BIOPROJECT: PRJEB21760

bioproject<-"PRJEB21760"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$treatment = relevel(as.factor(allsamples$treatment), ref="control")

#split by genotype and time, treat sex as blocking variable
for (i in c(2, 6, 10)) {
  for (j in c("resistant", "susceptible")) {
    bioproject2 <- paste0(bioproject, "-time", i, "-", j, collapse = "")
    subsample <- allsamples %>% filter(time == i, genotype == j)
    so<-load_sleuth_obj(subsample)
    export_pca(so, "treatment", bioproject2)
    export_pca(so, "sex", bioproject2)
    so<-sleuth_fit(so, ~treatment+sex, "treatment")
    export_results(so, "treatmentNDV", "treatment", bioproject2, "NDV")
  }
}

################

##BIOPROJECT: PRJNA143425

bioproject<-"PRJNA143425"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$timepoint = as.factor(allsamples$timepoint)
allsamples$timepoint = relevel(allsamples$timepoint, ref="0")
allsamples$treatment = as.factor(allsamples$treatment)
allsamples$treatment = relevel(allsamples$treatment, ref="control")
so<-load_sleuth_obj(allsamples)

#pca
export_pca(so, "treatment", bioproject)
export_pca(so, "timepoint", bioproject)

#analysis
#treatment plas vs control
so<-sleuth_fit(so, ~treatment, "treatment")
export_results(so, "treatmentMycoplasma_agalactiae", "treatment", bioproject, "Mycoplasma_agalactiae")

#timepoints vs 0
so<-sleuth_fit(so, ~timepoint, "timepoint")
export_results(so, "timepoint3", "timepoint", bioproject, "Mycoplasma_agalactiae")
export_results(so, "timepoint12", "timepoint", bioproject, "Mycoplasma_agalactiae")
export_results(so, "timepoint24", "timepoint", bioproject, "Mycoplasma_agalactiae")

###################

##BIOPROJECT: PRJNA174747

bioproject<-"PRJNA174747"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$treatment = relevel(as.factor(allsamples$treatment), ref="control")

so<-load_sleuth_obj(allsamples)
export_pca(so, "treatment", bioproject)
so<-sleuth_fit(so, ~treatment, "treatment")
export_results(so, "treatmentWNV", "treatment", bioproject, "WNV")

###################

##BIOPROJECT: PRJNA227801

bioproject<-"PRJNA227801"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$treatment = relevel(as.factor(allsamples$treatment), ref="MOCK")

#split by genotype and time, treat sex as blocking variable
for (i in c(2, 4)) {
    bioproject2 <- paste0(bioproject, "-time", i, collapse = "")
    subsample <- allsamples %>% filter(timepoint == i)
    so<-load_sleuth_obj(subsample)
    export_pca(so, "treatment", bioproject2)
    export_pca(so, "genotype", bioproject2)
    so<-sleuth_fit(so, ~treatment+genotype, "treatment")
    export_results(so, "treatmentMA15", "treatment", bioproject2, "SARS-MA15")
}

###################

##BIOPROJECT: PRJNA254331

bioproject<-"PRJNA254331"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$phenotype = as.factor(allsamples$phenotype)
allsamples$phenotype = relevel(allsamples$phenotype, ref="none")
allsamples$treatment = as.factor(allsamples$treatment)
allsamples$treatment = relevel(allsamples$treatment, ref="Control")
so<-load_sleuth_obj(allsamples)

#pca
export_pca(so, "treatment", bioproject)
export_pca(so, "phenotype", bioproject)

#analysis
#treatment plas vs control
so<-sleuth_fit(so, ~treatment, "treatment")
export_results(so, "treatmentEcoli", "treatment", bioproject, "Ecoli")

#timepoints vs 0
so<-sleuth_fit(so, ~phenotype, "phenotype")
export_results(so, "phenotypehigh_dose", "phenotype", bioproject, "Ecoli")
export_results(so, "phenotypelow_dose", "phenotype", bioproject, "Ecoli")
export_results(so, "phenotypesepsis", "phenotype", bioproject, "Ecoli")

###################

##BIOPROJECT: PRJNA288323

bioproject<-"PRJNA288323"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$phenotype = relevel(as.factor(allsamples$genotype), ref="control")

#susceptible and resistant are actually phenotypes, compare each to control at separate time points

so<-load_sleuth_obj(allsamples)

#split by time
for (i in c(1, 5)) {
  bioproject2 <- paste0(bioproject, "-time", i, collapse = "")
  subsample <- allsamples %>% filter(time == i)
  so<-load_sleuth_obj(subsample)
  export_pca(so, "treatment", bioproject2)
  export_pca(so, "phenotype", bioproject2)
  so<-sleuth_fit(so, ~phenotype, "phenotype")
  export_results(so, "phenotypesusceptible", "phenotype", bioproject2, "Ecoli")
  export_results(so, "phenotyperesistant", "phenotype", bioproject2, "Ecoli")
}

###################

##BIOPROJECT: PRJNA296924

bioproject<-"PRJNA296924"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$treatment = relevel(as.factor(allsamples$treatment), ref="control")

#split by genotype and time, treat sex as blocking variable
for (i in c(7, 28)) {
  for (j in c("L10L", "L10H")) {
    bioproject2 <- paste0(bioproject, "-time", i, "-", j, collapse = "")
    subsample <- allsamples %>% filter(time == i, genotype == j)
    so<-load_sleuth_obj(subsample)
    export_pca(so, "treatment", bioproject2)
    so<-sleuth_fit(so, ~treatment, "treatment")
    export_results(so, "treatmentIBV", "treatment", bioproject2, "IBV")
  }
}

################

##Bioproject PRJNA382632

bioproject<-"PRJNA382632"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$treatment = relevel(as.factor(allsamples$treatment), ref="control")

#split into timepoint, genotype as co-factor

for (i in c(3, 6, 12, 18)) {
    bioproject2 <- paste0(bioproject, "-time", i, collapse = "")
    subsample <- allsamples %>% filter(timepoint == i)
    so<-load_sleuth_obj(subsample)
    export_pca(so, "treatment", bioproject2)
    export_pca(so, "genotype", bioproject2)
    so<-sleuth_fit(so, ~treatment+genotype, "treatment")
    export_results(so, "treatmentH3N2", "treatment", bioproject2, "H3N2")
    export_results(so, "treatmentH5N1", "treatment", bioproject2, "H5N1")
    export_results(so, "treatmentH1N1", "treatment", bioproject2, "H1N1")
}

################



##Bioproject PRJNA394119

bioproject<-"PRJNA394119"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$treatment = relevel(as.factor(allsamples$treatment), ref="Control")

#split into timepoint, genotype as co-factor

for (i in c(1, 3, 5, 7)) {
  bioproject2 <- paste0(bioproject, "-time", i, collapse = "")
  subsample <- allsamples %>% filter(time == i)
  so<-load_sleuth_obj(subsample)
  export_pca(so, "treatment", bioproject2)
  so<-sleuth_fit(so, ~treatment, "treatment")
  export_results(so, "treatmentMgallisepticum", "treatment", bioproject2, "Mgallisepticum")
}

###################

##BIOPROJECT: PRJNA395845

bioproject<-"PRJNA395845"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$treatment = relevel(as.factor(allsamples$treatment), ref="control")

so<-load_sleuth_obj(allsamples)

#pca
export_pca(so, "treatment", bioproject)
export_pca(so, "delivery", bioproject)
export_pca(so, "sex", bioproject)

#analysis
#treatment plas vs control
so<-sleuth_fit(so, ~treatment+delivery+sex, "treatment")
export_results(so, "treatmentEcoli", "treatment", bioproject, "Ecoli")

################

##BIOPROJECT: PRJNA285798

bioproject<-"PRJNA285798"
allsamples<-load_samples_for_sleuth(bioproject)

#use treatment but restrict to eliminate healthy_children_with_bacterial_colonization
allsamples <- allsamples %>% filter(phenotype != "healthy_children_with_bacterial_colonization")
allsamples$treatment = relevel(as.factor(allsamples$treatment), ref="Negative")

so<-load_sleuth_obj(allsamples)

#pca
export_pca(so, "treatment", bioproject)
export_pca(so, "severity", bioproject)
export_pca(so, "age", bioproject)

#analysis
#treatment plas vs control
so<-sleuth_fit(so, ~treatment+age, "treatment")
export_results(so, "treatmentDAEC", "treatment", bioproject, "Ecoli")
export_results(so, "treatmentEAEC", "treatment", bioproject, "Ecoli")
export_results(so, "treatmentEPEC", "treatment", bioproject, "Ecoli")

################


##BIOPROJECT: PRJNA289121

bioproject<-"PRJNA289121"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$treatment = relevel(as.factor(allsamples$treatment), ref="mock")

#split into tissue, pool time points

for (j in c("lung", "brain")) {
    bioproject2 <- paste0(bioproject, "-", j, collapse = "")
    subsample <- allsamples %>% filter(tissue == j)
    so<-load_sleuth_obj(subsample)
    export_pca(so, "treatment", bioproject2)
    export_pca(so, "timepoint", bioproject2)    
    so<-sleuth_fit(so, ~treatment, "treatment")
    export_results(so, "treatmentHeV", "treatment", bioproject2, "HeV")
    export_results(so, "treatmentNiVB", "treatment", bioproject2, "NiVB")
}

#################

##BioProject PRJNA305099

bioproject<-"PRJNA305099"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$timepoint = as.factor(allsamples$timepoint)
allsamples$timepoint = relevel(allsamples$timepoint, ref="0")
allsamples$treatment = as.factor(allsamples$treatment)
allsamples$treatment = relevel(allsamples$treatment, ref="control")
so<-load_sleuth_obj(allsamples)

#pca
export_pca(so, "treatment", bioproject)
export_pca(so, "timepoint", bioproject)

#analysis
#treatment plas vs control
so<-sleuth_fit(so, ~treatment, "treatment")
export_results(so, "treatmentinfluenza", "treatment", bioproject, "influenza")

#timepoints vs 0
so<-sleuth_fit(so, ~timepoint, "timepoint")
export_results(so, "timepoint8", "timepoint", bioproject, "influenza")
export_results(so, "timepoint24", "timepoint", bioproject, "influenza")

################


##BioProject PRJNA314450

bioproject<-"PRJNA314450"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$timepoint = as.factor(allsamples$timepoint)
allsamples$timepoint = relevel(allsamples$timepoint, ref="0")
allsamples$treatment = as.factor(allsamples$treatment)
allsamples$treatment = relevel(allsamples$treatment, ref="Control")
so<-load_sleuth_obj(allsamples)

#pca
export_pca(so, "treatment", bioproject)
export_pca(so, "timepoint", bioproject)

#analysis
#treatment plas vs control
so<-sleuth_fit(so, ~treatment, "treatment")
export_results(so, "treatmentPlasmodium_berghei", "treatment", bioproject, "Plasmodium_berghei")

#timepoints vs 0
so<-sleuth_fit(so, ~timepoint, "timepoint")
export_results(so, "timepoint24", "timepoint", bioproject, "Plasmodium_berghei")
export_results(so, "timepoint48", "timepoint", bioproject, "Plasmodium_berghei")

##############


##BioProject PRJNA352507

bioproject<-"PRJNA352507"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$time = as.factor(allsamples$time)
allsamples$time = relevel(allsamples$time, ref="0")
allsamples$treatment = as.factor(allsamples$treatment)
allsamples$treatment = relevel(allsamples$treatment, ref="control")
so<-load_sleuth_obj(allsamples)

#pca
export_pca(so, "treatment", bioproject)
export_pca(so, "time", bioproject)

#analysis
#treatment plas vs control
so<-sleuth_fit(so, ~treatment, "treatment")
export_results(so, "treatmentWNV", "treatment", bioproject, "WNV")

#timepoints vs 0
so<-sleuth_fit(so, ~time, "time")
export_results(so, "time2", "time", bioproject, "WNV")
export_results(so, "time4", "time", bioproject, "WNV")


##############

##Bioproject PRJNA385346

bioproject<-"PRJNA385346"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$treatment = relevel(as.factor(allsamples$treatment), ref="control")

#pool all timepoints because of odd experimental design (controls are day 4)

so<-load_sleuth_obj(allsamples)
export_pca(so, "treatment", bioproject)
export_pca(so, "timepoint", bioproject)

so<-sleuth_fit(so, ~treatment, "treatment")
export_results(so, "treatmentH3N2", "treatment", bioproject, "H3N2")
export_results(so, "treatmentH5N1", "treatment", bioproject, "H5N1")
export_results(so, "treatmentH1N1", "treatment", bioproject, "H1N1")

###################


##BIOPROJECT: PRJNA257687

bioproject<-"PRJNA257687"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$time = relevel(as.factor(allsamples$time), ref="0")
#this one is a mess, just going to make new columns
#first is comparison of just bird2 - bird4, 21+31 vs 0
#second is comparison of all individuals infected vs not (pooling time 0 + control)
allsamples <- allsamples %>% 
  mutate(treatment2 = ifelse(treatment == "Avian_malaria" & time == "0", "Control", treatment)) %>%
  mutate(treatment3 = ifelse(treatment == "Avian_malaria" & time == "0", "Naive", treatment))
allsamples$treatment2 = relevel(as.factor(allsamples$treatment2), ref="Control")
allsamples$treatment3 = relevel(as.factor(allsamples$treatment3), ref="Naive")

so<-load_sleuth_obj(allsamples)

#pca
export_pca(so, "treatment", bioproject)
export_pca(so, "treatment2", bioproject)
export_pca(so, "treatment3", bioproject)
export_pca(so, "time", bioproject)

#analysis
#treatment plas vs control
so<-sleuth_fit(so, ~treatment2, "treatment2")
so<-sleuth_fit(so, ~treatment3, "treatment3")
export_results(so, "treatment2Avian_malaria", "treatment2", bioproject, "Avian_malaria")
export_results(so, "treatment3Avian_malaria", "treatment3", bioproject, "Avian_malaria")

###################


##BIOPROJECT: PRJNA279199

bioproject<-"PRJNA279199"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$treatment = relevel(as.factor(allsamples$treatment), ref="PRE")

so<-load_sleuth_obj(allsamples)

#pca
export_pca(so, "treatment", bioproject)
export_pca(so, "location", bioproject)
export_pca(so, "sex", bioproject)

#analysis
#treatment plas vs control
so<-sleuth_fit(so, ~treatment+location+sex, "treatment")
export_results(so, "treatmentDX", "treatment", bioproject, "Plasmodium_vivax")

###################


##BIOPROJECT: PRJNA279487

bioproject<-"PRJNA279487"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$phenotype = relevel(as.factor(allsamples$genotype), ref="control")

#susceptible and resistant are actually phenotypes, compare each to control at separate time points

so<-load_sleuth_obj(allsamples)

#split by time
for (i in c(1, 5)) {
  bioproject2 <- paste0(bioproject, "-time", i, collapse = "")
  subsample <- allsamples %>% filter(time == i)
  so<-load_sleuth_obj(subsample)
  export_pca(so, "treatment", bioproject2)
  export_pca(so, "phenotype", bioproject2)
  so<-sleuth_fit(so, ~phenotype, "phenotype")
  export_results(so, "phenotypesusceptible", "phenotype", bioproject2, "Ecoli")
  export_results(so, "phenotyperesistant", "phenotype", bioproject2, "Ecoli")
}


###################


##BIOPROJECT: PRJNA284293

bioproject<-"PRJNA284293"
allsamples<-load_samples_for_sleuth(bioproject)
allsamples$phenotype = relevel(as.factor(allsamples$genotype), ref="control")

#susceptible and resistant are actually phenotypes, compare each to control at separate time points

so<-load_sleuth_obj(allsamples)

#split by time
for (i in c(1, 5)) {
  bioproject2 <- paste0(bioproject, "-time", i, collapse = "")
  subsample <- allsamples %>% filter(time == i)
  so<-load_sleuth_obj(subsample)
  export_pca(so, "treatment", bioproject2)
  export_pca(so, "phenotype", bioproject2)
  so<-sleuth_fit(so, ~phenotype, "phenotype")
  export_results(so, "phenotypesusceptible", "phenotype", bioproject2, "Ecoli")
  export_results(so, "phenotyperesistant", "phenotype", bioproject2, "Ecoli")
}
