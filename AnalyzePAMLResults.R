setwd("~/Dropbox/BirdImmuneGeneEvolution")
require(ggplot2)

paml_res <- read.table("~/Dropbox/BirdImmuneGeneEvolution/PAML/paml_res_all_hogs.txt",header=T)

hog_ids <- read.table("~/Dropbox/BirdImmuneGeneEvolution/new_hog_list.txt")
colnames(hog_ids) <- c("HOG2_HogID","NCBI_ID","entrezgene","sp")

immune <- read.csv("~/Dropbox/BirdImmuneGeneEvolution/CuratedAvianImmuneGenes_NCBINums.csv",header=T,row.names=1)

innateDB <- read.csv("~/Dropbox/BirdImmuneGeneEvolution/InnateImmuneDB_NCBINums.csv")

#Separate out HOG IDs, extract only galGal sequences
hog_ids$hog <- apply(hog_ids,1,function(x) {unlist(strsplit(as.character(x[1]),"_"))})[2,]
hog_ids_gg <- hog_ids[grep("galGal",hog_ids[,4]),]

#Merge Hog IDs and NCBI gene IDs
paml_res_ncbi <- merge(hog_ids_gg,paml_res,by="hog")

paml_immfunc <- merge(paml_res_ncbi,immune,by="entrezgene")

funcplot <- ggplot(data=paml_immfunc,aes(x=Class,y=kappa,fill=Class)) + geom_boxplot() + geom_jitter()

funcplot + facet_grid(. ~ model_num) + theme_bw() + theme(axis.text.x=element_text(angle = -45, hjust = 0))

funcplot <- ggplot(data=paml_immfunc,aes(x=Class,y=treelen,fill=Class)) + geom_boxplot() + geom_jitter()

funcplot + facet_grid(. ~ model_num) + theme_bw() + theme(axis.text.x=element_text(angle = -45, hjust = 0))

#Function to conduct likelihood ratio tests between m1 and m2 models, and m7 and m8 models when given a hogid.
lrt_hog <- function(hogid,results){
	res_hog <- results[results$hog==hogid,]
	chisq_m1m2 <- 2*(res_hog[res_hog$model_num==2,"lnl"]-res_hog[res_hog$model_num==1,"lnl"])
	chisq_m7m8 <- 2*(res_hog[res_hog$model_num==8,"lnl"]-res_hog[res_hog$model_num==7,"lnl"])
	pval_m1m2 <- 1-pchisq(chisq_m1m2,2)
	pval_m7m8 <- 1-pchisq(chisq_m7m8,2)
	
	return(c(chisq_m1m2,pval_m1m2,chisq_m7m8,pval_m7m8))
}


hogs_with_ids <- unique(paml_res_ncbi$hog)
paml_pval_allgenes <- matrix(nrow=length(hogs_with_ids),ncol=4)

for (i in 1:length(unique(paml_res_ncbi))){
	paml_pval_allgenes[i,] <- lrt_hog(hogs_with_ids[i],paml_res)
}


set2 <- read.table("~/AvianImmuneRes/PAML/all_hogs_rewrite")
set2 <- read.table("~/AvianImmuneRes/PAML/test_hogs_100_2")
set2 <- read.table("~/AvianImmuneRes/PAML/split_hogs_3")
set2 <- read.table("~/AvianImmuneRes/PAML/split_hogs_2")
res2 <- read.table("~/AvianImmuneRes/PAML/paml_res_all_hogs_rewrite.txt",header=T)
res2 <- read.table("~/AvianImmuneRes/PAML/paml_res_test100.txt",header=T)
res2 <- read.table("~/AvianImmuneRes/PAML/paml_res_hogs3.txt",header=T)
res2 <- read.table("~/AvianImmuneRes/PAML/paml_res_hogs2.txt",header=T)
hos2 <- unique(res2$hog)

setdiff(as.character(set2[,1]),as.character(hos2))