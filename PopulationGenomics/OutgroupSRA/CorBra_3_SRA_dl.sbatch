#!/bin/bash
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-12:00 # Runtime in days-hours:minutes
#SBATCH --mem 16000 # Memory in MB
#SBATCH -J CorBra_3_dl # job name
#SBATCH -o CorBra_3_dl.log # File to which standard out will be written
#SBATCH -e CorBra_3_dl.err # File to which standard err will be written

wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR128/SRR1286117/SRR1286117.sra

wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR128/SRR1286118/SRR1286118.sra

wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR128/SRR1286119/SRR1286119.sra

wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR128/SRR1286121/SRR1286121.sra

wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR128/SRR1286122/SRR1286122.sra