# avian-immunity

A series of scripts for studying the evolution of immune genes across birds.

The first set of scripts extracts tree files (paml_tree_extract.py) from preivous PAML runs of the m0 model for use when fitting 1,2,7 and 8 site models, and then submits jobs to a cluster using slurm (paml_launcher_el.sh, run_paml_sites.sh, and submit_paml_jobs.sh). Some of these scripts are modifications from the ratite-genomics repository.

The second set of scripts parses the PAML output (site_parser.py).
