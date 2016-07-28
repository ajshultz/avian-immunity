setwd("~/Dropbox/BirdImmuneGeneEvolution")
require(biomaRt)

curGenes <- read.csv("~/Dropbox/BirdImmuneGeneEvolution/CuratedAvianImmuneGenes.csv")
innateGenes <- read.csv("~/Dropbox/BirdImmuneGeneEvolution/InnateImmuneDB.csv")


#Get the chicken gene mart
ensembl = useMart("ensembl")
ensembl = useDataset("ggallus_gene_ensembl",mart=ensembl)

attributes = listAttributes(ensembl)

#Note that a few IDs have multiple Entrez Gene IDs
convert_curGenes <- getBM(attributes=c("ensembl_gene_id","entrezgene","hgnc_id","hgnc_symbol"),filters="ensembl_gene_id",mart=ensembl,values=curGenes$Ens.ID)
convert_innateGenes <- getBM(attributes=c("ensembl_gene_id","entrezgene","hgnc_id","hgnc_symbol"),filters="ensembl_gene_id",mart=ensembl,values=innateGenes$ensembl_gene_id)

#Merge into a single dataframe
colnames(curGenes)[3] <-"ensembl_gene_id" 
curGenes_conv <- merge(curGenes,convert_curGenes)
innateGenes_conv <- merge(innateGenes,convert_innateGenes)

write.csv(curGenes_conv,file="~/Dropbox/BirdImmuneGeneEvolution/CuratedAvianImmuneGenes_NCBINums.csv")
write.csv(innateGenes_conv,file="~/Dropbox/BirdImmuneGeneEvolution/InnateImmuneDB_NCBINums.csv")

