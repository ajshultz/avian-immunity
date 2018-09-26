# avian-immunity

A series of scripts for studying the evolution of immune genes across birds and comparing these results to those in mammals from Enard et al. 2016.

**Comparative Genomics:**

1. PAML_Run: Extracts tree files (paml_tree_extract.py) from previous PAML runs of the m0 model for use when fitting 1,2,7 and 8 site models, and then submits jobs to the Harvard computing cluster using slurm (paml_launcher_el.sh, run_paml_sites.sh, and submit_paml_jobs.sh). Some of these scripts are modifications from the ratite-genomics repository.
	
  -Note that HOGs had to be split into multiple files to allow for max # jobs on cluster
		
  -codeml.sites.ctl file was sometimes modified to only include a specific model in cases where analyses ran too long, or when omega was fixed at 1 for 2a and 8a models.

2. PAML_Parse: Parses the PAML output (site_parser.py), and requires the branchsite_parser.py. Also parses HyPhy Busted output (busted_parser.py), aBS-REL output (overall = parse_hyphy_bs-rel.sh; to capture branch-specific p-values = parse_bsrel_hog_branch_pvals.py). Finally, extract_alignment_length.py caluclates alignment lengths from directory structure created with PAML_Run.
	
3. PAML_Analyze: 
  -bed_extract_ensembl_NCBIGeneID.py: Reads in a bed file with both ensembl IDs and NCBI gene IDs merged, and outputs a translation table of IDs
  -01_PAML_HyPhy_Res_DataPrep.R - Process all PAML and HyPhy output into analysis-ready formats.
  -02_GeneID_Annotation.R - Annotate all HOGs with NCBI and Ensembl gene IDs. 
  -03_Basic_Stats.R - General dataset statistic calculations.
  -04_Pathway_Enrichment.R - Bird KEGG pathway enrichment analyses.
  -05_Mammal_Bird_Gene_Set_Comparisons.R - Comparisons of bird and mammal selection scan results.
  -06_BSREL_species_clustering.R - Clustering and analysis of bird aBS-REL results.
  -07_Transciptomics.R - ranscriptome analyses for both birds and mammals.
		
**Transcriptomics:**

-RNA-seq_datasets_birds.csv = all details on bird RNAseq datasets
-RNA-seq_datasets_mammals.csv = all details on mammal RNAseq datasets
-01_dl_prep_indexes.sh - download all indexes used and process them with Kallisto
-02_download_process.py - Takes RNA-seq dataset files, and automates SRA downloading and Kallisto processing with the SLURM cluster architecture
-prep_data_for_sleuth.md - Workflow description for prepping files and loading into sleuth for differential expression analyses
-sleuth directory - all sleuth prep files for each bioproject, and R script to run them 
