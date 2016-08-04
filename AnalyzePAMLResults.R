setwd("~/Dropbox/BirdImmuneGeneEvolution")
require(ggplot2)

paml_res <- read.table("~/AvianImmuneRes/PAML/paml_res_allhogs_rewrite.txt",header=T)

hog_ids <- read.table("~/Dropbox/BirdImmuneGeneEvolution/new_hog_list.txt")
colnames(hog_ids) <- c("HOG2_HogID","NCBI_ID","entrezgene","sp")

immune <- read.csv("~/Dropbox/BirdImmuneGeneEvolution/CuratedAvianImmuneGenes_NCBINums.csv",header=T,row.names=1)
#214 genes have NCBI gene IDs/251

innateDB <- read.csv("~/Dropbox/BirdImmuneGeneEvolution/InnateImmuneDB_NCBINums.csv")
#838 genes have NCBI gene IDs/1050

#Separate out HOG IDs, extract only galGal sequences
hog_ids$hog <- apply(hog_ids,1,function(x) {unlist(strsplit(as.character(x[1]),"_"))})[2,]
hog_ids_gg <- hog_ids[grep("galGal",hog_ids[,4]),]
#17157 have NCBI gene IDs


#Merge Hog IDs and NCBI gene IDs
paml_res_ncbi <- merge(hog_ids_gg,paml_res,by="hog")
#45026 rows from 10906 hogs

#which hogs are not present in paml results?
setdiff(hog_ids_gg$hog,paml_res_ncbi$hog)
setdiff(paml_res_ncbi$hog,hog_ids_gg$hog)
#all PAML results have ncbi gene Ids

paml_immfunc <- merge(paml_res_ncbi,immune,by="entrezgene")
#497 matches from only 122 HOGs

funcplot <- ggplot(data=paml_immfunc,aes(x=Class,y=kappa,fill=Class)) + geom_boxplot() + geom_jitter()

funcplot + facet_grid(. ~ model_num) + theme_bw() + theme(axis.text.x=element_text(angle = -45, hjust = 0))

funcplot <- ggplot(data=paml_immfunc,aes(x=Class,y=treelen,fill=Class)) + geom_boxplot() + geom_jitter()

funcplot + facet_grid(. ~ model_num) + theme_bw() + theme(axis.text.x=element_text(angle = -45, hjust = 0))

#Function to conduct likelihood ratio tests between m1 and m2 models, and m7 and m8 models when given a hogid.
lrt_hog <- function(hogid,results){
	res_hog <- results[results$hog==hogid,]
	chisq_m1m2 <- chisq_m7m8 <- pval_m1m2 <- pval_m7m8 <- NA
	
	if (1 %in% res_hog$model_num & 2 %in% res_hog$model_num){
		chisq_m1m2 <- 2*(res_hog[res_hog$model_num==2,"lnl"]-res_hog[res_hog$model_num==1,"lnl"])
		pval_m1m2 <- 1-pchisq(chisq_m1m2,2)
	}
	
	if (7 %in% res_hog$model_num & 8 %in% res_hog$model_num){
		chisq_m7m8 <- 2*(res_hog[res_hog$model_num==8,"lnl"]-res_hog[res_hog$model_num==7,"lnl"])
		pval_m7m8 <- 1-pchisq(chisq_m7m8,2)
	}
	
	
	return(c(chisq_m1m2,pval_m1m2,chisq_m7m8,pval_m7m8))
}

#Calculate p-values for all hogs with NCBI IDs between neutral and positive selection models (including FDR p-values)
hogs_with_ids <- unique(paml_res_ncbi$hog)
paml_pval_allgenes <- matrix(nrow=length(hogs_with_ids),ncol=4)

for (i in 1:length(hogs_with_ids)){
	paml_pval_allgenes[i,] <- lrt_hog(hogs_with_ids[i],paml_res)
}
paml_pval_allgenes <- data.frame(paml_pval_allgenes)
colnames(paml_pval_allgenes) <- c("ChiSq_m1m2","PVal_m1m2","ChiSq_m7m8","PVal_m7m8")
rownames(paml_pval_allgenes) <- hogs_with_ids

paml_pval_allgenes$FDRPval_m1m2 <- p.adjust(paml_pval_allgenes$PVal_m1m2,method="BH")
paml_pval_allgenes$FDRPval_m7m8 <- p.adjust(paml_pval_allgenes$PVal_m7m8,method="BH")

paml_pval_allgenes$hog <- rownames(paml_pval_allgenes)

paml_pval_allgenes_ncbi <- merge(hog_ids_gg,paml_pval_allgenes,by="hog")
paml_pval_allgenes_ncbi$entrezgene <- as.character(paml_pval_allgenes_ncbi$entrezgene)

write.csv(paml_pval_allgenes_ncbi,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/allgenePAML_Pvals.csv")

#Write file with all entrezgenes for all sig. m1m2 reactions.
allgenem1m2 <- paml_pval_allgenes_ncbi[paml_pval_allgenes_ncbi $FDRPval_m1m2<0.05,"entrezgene"]
allgenem1m2 <- allgenem1m2[!is.na(allgenem1m2)]
write.csv(allgenem1m2,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/EntrezGeneID_Allm1m2SigGenes.csv",row.names=F,quote=F)


immune$match <- immune$entrezgene %in% paml_pval_allgenes_ncbi$entrezgene
#Which immune genes without PAML results?
rownames(hog_ids_gg) <- hog_ids_gg$entrezgene
hog_ids_gg[as.character(immune$entrezgene[!immune$match][!is.na(immune$entrezgene[!immune$match])]),]
#Looking through the PAML folders, it appears that the discrepancy is due to no PAML models actually being available for those HOGS (must not have passed filters). Some of the entrezgenes match to the same hog (e.g. 882 shows up a number of times.)

paml_pval_immune <- merge(paml_pval_allgenes_ncbi,immune,by="entrezgene",incomparables=NA)
#126 genes have an annotated immune function from Julia's work, and PAML results. There are 42 effector genes, 46 receptor genes, 35 signaling genes, and 3 effector/signaling

paml_pval_innatedb <- merge(paml_pval_allgenes_ncbi,innateDB,by="entrezgene",incomparables=NA)

write.csv(paml_pval_immune,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/Paml_Pvals_CuratedAvianImmuneGenes.csv")
write.csv(paml_pval_innatedb,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/Paml_Pvals_InnateImmuneDB.csv")

#What proportion of immune classes are positively selected?
immunecounts <- as.data.frame(table(paml_pval_immune$Class))
rownames(immunecounts) <- immunecounts$Var1
immunecounts$Var1 <- NULL
colnames(immunecounts) <- c("NumGenes")

for (i in 1:nrow(immunecounts)){
	immunecounts$m1m2select[i] <- nrow(paml_pval_immune[paml_pval_immune$FDRPval_m1m2<0.05 & paml_pval_immune$Class==rownames(immunecounts)[i],])
	immunecounts$m7m8select[i] <- nrow(paml_pval_immune[paml_pval_immune$FDRPval_m7m8<0.05 & paml_pval_immune$Class==rownames(immunecounts)[i],])
}

#Add in results/expectations from all genes, all immune genes Julia curated, and innate immuneDB genes
immunecounts["AllGenes",] <- c(nrow(paml_pval_allgenes),nrow(paml_pval_allgenes[paml_pval_allgenes$FDRPval_m1m2<0.05,]),nrow(paml_pval_allgenes[paml_pval_allgenes$FDRPval_m7m8<0.05,]))

immunecounts["ImmuneGenes",] <- c(nrow(paml_pval_immune),nrow(paml_pval_immune[paml_pval_immune$FDRPval_m1m2<0.05,]),nrow(paml_pval_immune[paml_pval_immune$FDRPval_m7m8<0.05,]))

immunecounts["InnateDBGenes",] <- c(nrow(paml_pval_innatedb),nrow(paml_pval_innatedb[paml_pval_innatedb$FDRPval_m1m2<0.05,]),nrow(paml_pval_innatedb[paml_pval_innatedb$FDRPval_m7m8<0.05,]))

immunecounts$prop_m1m2select <- round(immunecounts$m1m2select/immunecounts$NumGenes,digits=3)
immunecounts$prop_m7m8select <- round(immunecounts$m7m8select/immunecounts$NumGenes,digits=3)

write.csv(immunecounts,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/InitialSignificanceResults.csv")

immunecounts$m1m2notselect <- immunecounts$NumGenes-immunecounts$m1m2select
immunecounts$m7m8notselect <- immunecounts$NumGenes-immunecounts$m7m8select

#Fisher's Exact test for signficance of different immune classes vs. all immune genes with the m1 vs. m2 model
#Effector vs. immune genes
fisher.test(immunecounts[c("Effector","ImmuneGenes"),c("m1m2select","m1m2notselect")])
#Receptor vs. immune genes
fisher.test(immunecounts[c("Receptor","ImmuneGenes"),c("m1m2select","m1m2notselect")])
#Signaling vs. immune genes
fisher.test(immunecounts[c("Signaling","ImmuneGenes"),c("m1m2select","m1m2notselect")])
#None are signficant

#Fisher's Exact test for signficance of different immune classes vs. all immune genes with the m7 vs. m8 model
#Effector vs. immune genes
fisher.test(immunecounts[c("Effector","ImmuneGenes"),c("m7m8select","m7m8notselect")])
#Receptor vs. immune genes
fisher.test(immunecounts[c("Receptor","ImmuneGenes"),c("m7m8select","m7m8notselect")])
#Signaling vs. immune genes
fisher.test(immunecounts[c("Signaling","ImmuneGenes"),c("m7m8select","m7m8notselect")])

barplot(immunecounts$prop_m1m2select,names.arg=rownames(immunecounts))
barplot(immunecounts$prop_m7m8select,names.arg=rownames(immunecounts))
