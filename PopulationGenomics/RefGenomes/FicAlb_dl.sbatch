#!/bin/bash
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-2:00 # Runtime in days-hours:minutes
#SBATCH --mem 1000 # Memory in MB
#SBATCH -J FicAlb_dl # job name
#SBATCH -o FicAlb_dl.out # File to which standard out will be written
#SBATCH -e FicAlb_dl.err # File to which standard err will be written

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000247815.1_FicAlb1.5/GCF_000247815.1_FicAlb1.5_genomic.fna.gz

cp GCF_000247815.1_FicAlb1.5_genomic.fna.gz FicAlb.fna.gz