#!/bin/bash

wget ftp://ftp.ensembl.org/pub/release-92/fasta/taeniopygia_guttata/cdna/README

#Capra hiracus
wget ftp://ftp.ensembl.org/pub/release-92/fasta/capra_hircus/cdna/Capra_hircus.ARS1.cdna.all.fa.gz
mv Capra_hircus.ARS1.cdna.all.fa.gz Capra_hircus.cdna.fa.gz
gunzip Capra_hircus.cdna.fa.gz

#Homo sapiens
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
mv Homo_sapiens.GRCh38.cdna.all.fa.gz Homo_sapiens.cdna.fa.gz
gunzip Homo_sapiens.cdna.fa.gz

#Mustela putoris furo
wget ftp://ftp.ensembl.org/pub/release-92/fasta/mustela_putorius_furo/cdna/Mustela_putorius_furo.MusPutFur1.0.cdna.all.fa.gz
mv Mustela_putorius_furo.MusPutFur1.0.cdna.all.fa.gz Mustela_putorius_furo.cdna.fa.gz
gunzip Mustela_putorius_furo.cdna.fa.gz

#Mus musculus
wget ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
mv Mus_musculus.GRCm38.cdna.all.fa.gz Mus_musculus.cdna.fa.gz
gunzip Mus_musculus.cdna.fa.gz

#Macaca fascicularis
wget ftp://ftp.ensembl.org/pub/release-92/fasta/macaca_fascicularis/cdna/Macaca_fascicularis.Macaca_fascicularis_5.0.cdna.all.fa.gz
mv Macaca_fascicularis.Macaca_fascicularis_5.0.cdna.all.fa.gz Macaca_fascicularis.cdna.fa.gz
gunzip Macaca_fascicularis.cdna.fa.gz

#Ovis aries
wget ftp://ftp.ensembl.org/pub/release-92/fasta/ovis_aries/cdna/Ovis_aries.Oar_v3.1.cdna.all.fa.gz
mv Ovis_aries.Oar_v3.1.cdna.all.fa.gz Ovis_aries.cdna.fa.gz
gunzip Ovis_aries.cdna.fa.gz

#Taeniopygia guttata
wget ftp://ftp.ensembl.org/pub/release-92/fasta/taeniopygia_guttata/cdna/Taeniopygia_guttata.taeGut3.2.4.cdna.all.fa.gz
mv Taeniopygia_guttata.taeGut3.2.4.cdna.all.fa.gz Taeniopygia_guttata.cdna.fa.gz
gunzip Taeniopygia_guttata.cdna.fa.gz


#Gallus gallus
wget ftp://ftp.ensembl.org/pub/release-92/fasta/gallus_gallus/cdna/Gallus_gallus.Gallus_gallus-5.0.cdna.all.fa.gz
mv Gallus_gallus.Gallus_gallus-5.0.cdna.all.fa.gz Gallus_gallus.cdna.fa.gz
gunzip Gallus_gallus.cdna.fa.gz


#Anas platyrhynchos
wget ftp://ftp.ensembl.org/pub/release-92/fasta/anas_platyrhynchos/cdna/Anas_platyrhynchos.BGI_duck_1.0.cdna.all.fa.gz
mv Anas_platyrhynchos.BGI_duck_1.0.cdna.all.fa.gz Anas_platyrhynchos.cdna.fa.gz
gunzip Anas_platyrhynchos.cdna.fa.gz

#Serinus canaria
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/534/875/GCF_000534875.1_SCA1/GCF_000534875.1_SCA1_rna.fna.gz
mv GCF_000534875.1_SCA1_rna.fna.gz Serinus_canaria.cdna.fa.gz
gunzip Serinus_canaria.cdna.fa.gz



#Index with Kallisto
module load hdf5/1.8.12-fasrc08
module load gcc/4.8.2-fasrc01 kallisto/0.43.1-fasrc01
screen -S kallisto_index

for SPECIES in Capra_hircus Homo_sapiens Mustela_putorius_furo Mus_musculus  Macaca_fascicularis Ovis_aries Taeniopygia_guttata Gallus_gallus Anas_platyrhynchos Serinus_canaria
do
kallisto index -i ${SPECIES} ${SPECIES}.cdna.fa 
done

