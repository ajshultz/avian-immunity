setwd("~/Dropbox/BirdImmuneGeneEvolution")
require(ggplot2)

paml_res <- read.table("~/AvianImmuneRes/PAML/paml_res_allhogs_rewrite.txt",header=T)

hog_ids <- read.table("~/Dropbox/BirdImmuneGeneEvolution/new_hog_list.txt")
colnames(hog_ids) <- c("HOG2_HogID","NCBI_ID","entrezgene","sp")

immune <- read.csv("~/Dropbox/BirdImmuneGeneEvolution/CuratedAvianImmuneGenes_NCBINums.csv",header=T,row.names=1)
#214 genes have NCBI gene IDs/251

#Integrate NCBI numbers from Biomart vs. galGal positions:
for (i in 1:nrow(immune)){
	if (is.na(immune[i,"entrezgene"])){
		immune[i,"entrezgene"] <- immune[i,"NCBI_entrez_galGal"]
	}
}
#Now all genes have NCBI gene IDs

innateDB <- read.csv("~/Dropbox/BirdImmuneGeneEvolution/InnateImmuneDB_NCBINums.csv")
#838 genes have NCBI gene IDs/1050

#Integrate NCBI numbers from Biomart vs. galGal positions:
for (i in 1:nrow(innateDB)){
	if (is.na(innateDB[i,"entrezgene"])){
		innateDB[i,"entrezgene"] <- innateDB[i,"NCBI_entrez_galGal"]
	}
}
nrow(innateDB[is.na(innateDB$entrezgene),])

#Now only 13 are missing NCBI gene IDs


#Separate out HOG IDs, extract only galGal sequences
hog_ids$hog <- apply(hog_ids,1,function(x) {unlist(strsplit(as.character(x[1]),"_"))})[2,]
hog_ids_gg <- hog_ids[grep("galGal",hog_ids[,4]),]
#17157 have NCBI gene IDs


#Merge Hog IDs and NCBI gene IDs
paml_res_ncbi <- merge(hog_ids_gg,paml_res,by="hog")
#45026 rows from 10906 hogs

#which hogs are not present in paml results?
#setdiff(hog_ids_gg$hog,paml_res_ncbi$hog)
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

#This function will take teh results for a given hog, and extract the omega values and proportion of sites that fall under that omega for the m2 model.
omega_m2 <- function(hogid,results){
	res_hog <- results[results$hog==hogid,]
	prop_m2 <- omega_m2 <- NA
	
	if (2 %in% res_hog$model_num){
		omega_m2_hog <- as.character(res_hog[res_hog$model_num==2,"omega"])
		omega_m2_hog <- strsplit(omega_m2_hog,split=",")[[1]]
		omega_m2_hog <- omega_m2_hog[length(omega_m2_hog)]
		omega_m2_hog <- strsplit(omega_m2_hog,split=":")[[1]]
		prop_m2 <- as.numeric(omega_m2_hog[1])
		omega_m2 <- as.numeric(omega_m2_hog[2])
		}
	
	return(c(prop_m2,omega_m2))
}

#This function will take teh results for a given hog, and extract the omega values and proportion of sites that fall under that omega for the m8 model.
omega_m8 <- function(hogid,results){
	res_hog <- results[results$hog==hogid,]
	prop_m8 <- omega_m8 <- NA
	
	if (8 %in% res_hog$model_num){
		omega_m8_hog <- as.character(res_hog[res_hog$model_num==8,"omega"])
		omega_m8_hog <- strsplit(omega_m8_hog,split=",")[[1]]
		omega_m8_hog <- omega_m8_hog[length(omega_m8_hog)]
		omega_m8_hog <- strsplit(omega_m8_hog,split=":")[[1]]
		prop_m8 <- as.numeric(omega_m8_hog[1])
		omega_m8 <- as.numeric(omega_m8_hog[2])
		}
	
	return(c(prop_m8,omega_m8))
}

#Calculate p-values for all hogs with NCBI IDs between neutral and positive selection models (including FDR p-values)
hogs_with_ids <- unique(paml_res_ncbi$hog)
paml_pval_allgenes <- matrix(nrow=length(hogs_with_ids),ncol=8)

for (i in 1:length(hogs_with_ids)){
	paml_pval_allgenes[i,1:4] <- lrt_hog(hogs_with_ids[i],paml_res)
	paml_pval_allgenes[i,5:6] <- try(omega_m2(hogs_with_ids[i],paml_res))
	paml_pval_allgenes[i,7:8] <- try(omega_m8(hogs_with_ids[i],paml_res))
}
paml_pval_allgenes <- data.frame(paml_pval_allgenes)
colnames(paml_pval_allgenes) <- c("ChiSq_m1m2","PVal_m1m2","ChiSq_m7m8","PVal_m7m8","Prop_m2","Omega_m2","Prop_m8","Omega_m8")
rownames(paml_pval_allgenes) <- hogs_with_ids

paml_pval_allgenes$FDRPval_m1m2 <- p.adjust(paml_pval_allgenes$PVal_m1m2,method="BH")
paml_pval_allgenes$FDRPval_m7m8 <- p.adjust(paml_pval_allgenes$PVal_m7m8,method="BH")

paml_pval_allgenes$hog <- rownames(paml_pval_allgenes)

#How many hogs didn't finish model m2 or m8?
nrow(paml_pval_allgenes[is.na(paml_pval_allgenes$Omega_m2),])
nrow(paml_pval_allgenes[is.na(paml_pval_allgenes$Omega_m8),])
#11 hogs don't have m2 models.
#1001 hogs don't have m8 models.

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
#hog_ids_gg[as.character(immune$entrezgene[!immune$match][!is.na(immune$entrezgene[!immune$match])]),]
#Looking through the PAML folders, it appears that the discrepancy is due to no PAML models actually being available for those HOGS (must not have passed filters). Some of the entrezgenes match to the same hog (e.g. 882 shows up a number of times.)

paml_pval_immune <- merge(paml_pval_allgenes_ncbi,immune,by="entrezgene",incomparables=NA)
#145 genes have an annotated immune function from Julia's work, and PAML results. There are 48 effector genes, 51 receptor genes, 43 signaling genes, and 3 effector/signaling

#length(unique(as.character(paml_pval_immune$Name)))
#There seem to be some multiples, need to filter those
paml_pval_immune[paml_pval_immune$Name %in% names(table(paml_pval_immune$Name))[table(paml_pval_immune$Name)>1],]

#Command to filter duplicates:
paml_pval_immune <- paml_pval_immune[!duplicated(paml_pval_immune$Name),]

paml_pval_innatedb <- merge(paml_pval_allgenes_ncbi,innateDB,by="entrezgene",incomparables=NA)

#Also remove duplicates from innateDB list.
paml_pval_innatedb <- paml_pval_innatedb[!duplicated(paml_pval_innatedb$InnateDBID),]
#694 results for the innatedb.

write.csv(paml_pval_immune,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/Paml_Pvals_CuratedAvianImmuneGenes.csv")
write.csv(paml_pval_innatedb,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/Paml_Pvals_InnateImmuneDB.csv")



#What proportion of immune classes are positively selected?
immunecounts <- as.data.frame(table(paml_pval_immune$Class))
rownames(immunecounts) <- immunecounts$Var1
immunecounts$Var1 <- NULL
colnames(immunecounts) <- c("NumGenes")

for (i in 1:nrow(immunecounts)){
	immunecounts$m1m2select[i] <- nrow(paml_pval_immune[paml_pval_immune$FDRPval_m1m2<0.05 & paml_pval_immune$Class==rownames(immunecounts)[i],])
	immunecounts$m1m2selectOmega[i] <- nrow(paml_pval_immune[paml_pval_immune$FDRPval_m1m2<0.05 & paml_pval_immune$Class==rownames(immunecounts)[i] & paml_pval_immune$Omega_m2>2 & paml_pval_immune$Prop_m2>0.05,])
	immunecounts$m7m8select[i] <- nrow(paml_pval_immune[paml_pval_immune$FDRPval_m7m8<0.05 & paml_pval_immune$Class==rownames(immunecounts)[i],])
	immunecounts$m7m8selectOmega[i] <- nrow(paml_pval_immune[paml_pval_immune$FDRPval_m7m8<0.05 & paml_pval_immune$Class==rownames(immunecounts)[i] & paml_pval_immune$Omega_m8>2 & paml_pval_immune$Prop_m8>0.05,])
	immunecounts$m2modelpresent[i] <- nrow(paml_pval_immune[!is.na(paml_pval_immune$Omega_m2) & paml_pval_immune$Class==rownames(immunecounts)[i],])
	immunecounts$m8modelpresent[i] <- nrow(paml_pval_immune[!is.na(paml_pval_immune$Omega_m8) & paml_pval_immune$Class==rownames(immunecounts)[i],])
}

#Add in results/expectations from all genes, all immune genes Julia curated, and innate immuneDB genes
immunecounts["AllGenes",] <- c(nrow(paml_pval_allgenes),nrow(paml_pval_allgenes[paml_pval_allgenes$FDRPval_m1m2<0.05,]),nrow(paml_pval_allgenes[paml_pval_allgenes$FDRPval_m1m2<0.05 & paml_pval_allgenes$Omega_m2>2 & paml_pval_allgenes$Prop_m2>0.05,]),nrow(paml_pval_allgenes[paml_pval_allgenes$FDRPval_m7m8<0.05,]),nrow(paml_pval_allgenes[paml_pval_allgenes$FDRPval_m7m8<0.05 & paml_pval_allgenes$Omega_m8>2 & paml_pval_allgenes$Prop_m8>0.05,]),nrow(paml_pval_allgenes[!is.na(paml_pval_allgenes$Omega_m2),]),nrow(paml_pval_allgenes[!is.na(paml_pval_allgenes$Omega_m8),]))

immunecounts["ImmuneGenes",] <- c(nrow(paml_pval_immune),nrow(paml_pval_immune[paml_pval_immune$FDRPval_m1m2<0.05,]),nrow(paml_pval_immune[paml_pval_immune$FDRPval_m1m2<0.05 & paml_pval_immune$Omega_m2>2 & paml_pval_immune$Prop_m2>0.05,]),nrow(paml_pval_immune[paml_pval_immune$FDRPval_m7m8<0.05,]),nrow(paml_pval_immune[paml_pval_immune$FDRPval_m7m8<0.05 & paml_pval_immune$Omega_m8>2 & paml_pval_immune$Prop_m8>0.05,]),nrow(paml_pval_immune[!is.na(paml_pval_immune$Omega_m2),]),nrow(paml_pval_immune[!is.na(paml_pval_immune$Omega_m8),]))

immunecounts["InnateDBGenes",] <- c(nrow(paml_pval_innatedb),nrow(paml_pval_innatedb[paml_pval_innatedb$FDRPval_m1m2<0.05,]),nrow(paml_pval_innatedb[paml_pval_innatedb$FDRPval_m1m2<0.05 & paml_pval_innatedb$Omega_m2>2 & paml_pval_innatedb$Prop_m2>0.05,]),nrow(paml_pval_innatedb[paml_pval_innatedb$FDRPval_m7m8<0.05,]),nrow(paml_pval_innatedb[paml_pval_innatedb$FDRPval_m7m8<0.05 & paml_pval_innatedb$Omega_m8>2 & paml_pval_innatedb$Prop_m8>0.05,]),nrow(paml_pval_innatedb[!is.na(paml_pval_innatedb$Omega_m2),]),nrow(paml_pval_innatedb[!is.na(paml_pval_innatedb$Omega_m8),]))

immunecounts$prop_m1m2select <- round(immunecounts$m1m2select/immunecounts$m2modelpresent,digits=3)
immunecounts$prop_m7m8select <- round(immunecounts$m7m8select/immunecounts$m8modelpresent,digits=3)
immunecounts$prop_m1m2selectOmega <- round(immunecounts$m1m2selectOmega/immunecounts$m2modelpresent,digits=3)
immunecounts$prop_m7m8selectOmega <- round(immunecounts$m7m8selectOmega/immunecounts$m8modelpresent,digits=3)


write.csv(immunecounts,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/InitialSignificanceResults.csv")

immunecounts$m1m2notselect <- immunecounts$m2modelpresent-immunecounts$m1m2select
immunecounts$m7m8notselect <- immunecounts$m8modelpresent-immunecounts$m7m8select
immunecounts$m1m2notselectOmega <- immunecounts$m2modelpresent-immunecounts$m1m2selectOmega
immunecounts$m7m8notselectOmega <- immunecounts$m8modelpresent-immunecounts$m7m8selectOmega



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


#Fisher's Exact test for signficance of different immune classes vs. all immune genes with the m1 vs. m2 model with omega > 2 and with at least 5% of sites
#Effector vs. immune genes
fisher.test(immunecounts[c("Effector","ImmuneGenes"),c("m1m2selectOmega","m1m2notselectOmega")])
#Receptor vs. immune genes
fisher.test(immunecounts[c("Receptor","ImmuneGenes"),c("m1m2selectOmega","m1m2notselectOmega")])
#Signaling vs. immune genes
fisher.test(immunecounts[c("Signaling","ImmuneGenes"),c("m1m2selectOmega","m1m2notselectOmega")])
#None are signficant

#Fisher's Exact test for signficance of different immune classes vs. all immune genes with the m7 vs. m8 model with omega > 2 and with at least 5% of sites
#Effector vs. immune genes
fisher.test(immunecounts[c("Effector","ImmuneGenes"),c("m7m8selectOmega","m7m8notselectOmega")])
#Receptor vs. immune genes
fisher.test(immunecounts[c("Receptor","ImmuneGenes"),c("m7m8selectOmega","m7m8notselectOmega")])
#Signaling vs. immune genes
fisher.test(immunecounts[c("Signaling","ImmuneGenes"),c("m7m8selectOmega","m7m8notselectOmega")])

barplot(immunecounts$prop_m1m2select,names.arg=rownames(immunecounts))
barplot(immunecounts$prop_m7m8select,names.arg=rownames(immunecounts))
