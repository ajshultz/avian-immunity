#!/bin/bash
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-12:00 # Runtime in days-hours:minutes
#SBATCH --mem 16000 # Memory in MB
#SBATCH -J CorFru_1_dl # job name
#SBATCH -o CorFru_1_dl.log # File to which standard out will be written
#SBATCH -e CorFru_1_dl.err # File to which standard err will be written

wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR136/SRR1368309/SRR1368309.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR136/SRR1368322/SRR1368322.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR136/SRR1368329/SRR1368329.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR136/SRR1368330/SRR1368330.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR136/SRR1368331/SRR1368331.sra