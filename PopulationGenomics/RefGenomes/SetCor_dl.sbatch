#!/bin/bash
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-2:00 # Runtime in days-hours:minutes
#SBATCH --mem 1000 # Memory in MB
#SBATCH -J SetCor_dl # job name
#SBATCH -o SetCor_dl.out # File to which standard out will be written
#SBATCH -e SetCor_dl.err # File to which standard err will be written


wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/746/935/GCA_001746935.1_mywagenomev1.1/GCA_001746935.1_mywagenomev1.1_genomic.fna.gz

cp GCA_001746935.1_mywagenomev1.1_genomic.fna.gz SetCor.fna.gz