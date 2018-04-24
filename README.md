# avian-immunity

A series of scripts for studying the evolution of immune genes across birds. This is a work in progress, so all scripts are subject to frequent updating.

**Comparative Genomics:**

1. PAML_Run: Extracts tree files (paml_tree_extract.py) from previous PAML runs of the m0 model for use when fitting 1,2,7 and 8 site models, and then submits jobs to the Harvard computing cluster using slurm (paml_launcher_el.sh, run_paml_sites.sh, and submit_paml_jobs.sh). Some of these scripts are modifications from the ratite-genomics repository.
	
  -Note that HOGs had to be split into multiple files to allow for max # jobs on cluster
		
  -codeml.sites.ctl file was sometimes modified to only include a specific model in cases where analyses ran too long, or when omega was fixed at 1 for 2a and 8a models.

2. PAML_Parse: Parses the PAML output (site_parser.py), and requires the branchsite_parser.py created by Tim Sackton, modified to use Python 2.
	
3. PAML_Analyze: 
  -bed_extract_ensembl_NCBIGeneID.py: Reads in a bed file with both ensembl IDs and NCBI gene IDs merged, and outputs a translation table of IDs
		
  -EnsemblNumToNCBIGeneIDs.R: R script to convert ensembl IDs from curated lists of immune genes to NCBI Gene IDs, both using those available in Ensembl and from those extracted from bed file (above)
		
  -hogsToRerun.R: Quick R script to take parsed PAML results and identify HOGs that do not have all PAML models run
		
**Transcriptomics:**

