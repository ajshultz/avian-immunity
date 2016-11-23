#!/bin/bash

#Indexes all downloaded reference genomes for pop genomic projects.

sbatch /n/regal/edwards_lab/hofi/AvianImmune/avian-immunity/PopulationGenomics/RefGenomes/make_bwa_index.sh CorCor
sbatch /n/regal/edwards_lab/hofi/AvianImmune/avian-immunity/PopulationGenomics/RefGenomes/make_bwa_index.sh FicAlb
sbatch /n/regal/edwards_lab/hofi/AvianImmune/avian-immunity/PopulationGenomics/RefGenomes/make_bwa_index.sh ParMaj
sbatch /n/regal/edwards_lab/hofi/AvianImmune/avian-immunity/PopulationGenomics/RefGenomes/make_bwa_index.sh SetCor