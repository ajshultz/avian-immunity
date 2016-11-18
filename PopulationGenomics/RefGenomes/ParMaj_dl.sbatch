#!/bin/bash
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-2:00 # Runtime in days-hours:minutes
#SBATCH --mem 1000 # Memory in MB
#SBATCH -J ParMaj_dl # job name
#SBATCH -o ParMaj_dl.out # File to which standard out will be written
#SBATCH -e ParMaj_dl.err # File to which standard err will be written


wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/522/545/GCF_001522545.1_Parus_major1.0.3/GCF_001522545.1_Parus_major1.0.3_genomic.fna.gz

cp GCF_001522545.1_Parus_major1.0.3_genomic.fna.gz ParMaj.fna.gz