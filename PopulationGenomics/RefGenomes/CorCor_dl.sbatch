#!/bin/bash
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-2:00 # Runtime in days-hours:minutes
#SBATCH --mem 1000 # Memory in MB
#SBATCH -J CorCor_dl # job name
#SBATCH -o CorCor_dl.out # File to which standard out will be written
#SBATCH -e CorCor_dl.err # File to which standard err will be written


wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/738/735/GCF_000738735.1_Hooded_Crow_genome/GCF_000738735.1_Hooded_Crow_genome_genomic.fna.gz

cp GCF_000738735.1_Hooded_Crow_genome_genomic.fna.gz CorCor.fna.gz