setwd("~/Dropbox/BirdImmuneGeneEvolution")
require(ggplot2)

paml_res <- read.table("~/Dropbox/BirdImmuneGeneEvolution/PAML/paml_res_allhogs_rewrite.txt",header=T)

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

save(paml_pval_allgenes,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/paml_pvals_allgenes.RDat")

#How many hogs didn't finish model m2 or m8?
nrow(paml_pval_allgenes[is.na(paml_pval_allgenes$Omega_m2),])
nrow(paml_pval_allgenes[is.na(paml_pval_allgenes$Omega_m8),])
#11 hogs don't have m2 models.
#1001 hogs don't have m8 models.

paml_pval_allgenes_ncbi <- merge(hog_ids_gg,paml_pval_allgenes,by="hog")
paml_pval_allgenes_ncbi$entrezgene <- as.character(paml_pval_allgenes_ncbi$entrezgene)

#write.csv(paml_pval_allgenes_ncbi,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/allgenePAML_Pvals.csv")

#Write file with all entrezgenes for all sig. m1m2 reactions.
allgenem1m2 <- paml_pval_allgenes_ncbi[paml_pval_allgenes_ncbi $FDRPval_m1m2<0.05,"entrezgene"]
allgenem1m2 <- allgenem1m2[!is.na(allgenem1m2)]
#write.csv(allgenem1m2,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/EntrezGeneID_Allm1m2SigGenes.csv",row.names=F,quote=F)


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

#write.csv(paml_pval_immune,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/Paml_Pvals_CuratedAvianImmuneGenes.csv")
#write.csv(paml_pval_innatedb,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/Paml_Pvals_InnateImmuneDB.csv")



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


#write.csv(immunecounts,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/InitialSignificanceResults.csv")

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




#Versus all genes?
#Fisher's Exact test for signficance of different immune classes vs. all immune genes with the m1 vs. m2 model
#Effector vs. immune genes
fisher.test(immunecounts[c("Effector","AllGenes"),c("m1m2select","m1m2notselect")])
#Receptor vs. immune genes
fisher.test(immunecounts[c("Receptor","AllGenes"),c("m1m2select","m1m2notselect")])
#Signaling vs. immune genes
fisher.test(immunecounts[c("Signaling","AllGenes"),c("m1m2select","m1m2notselect")])
#None are signficant

#Fisher's Exact test for signficance of different immune classes vs. all immune genes with the m7 vs. m8 model
#Effector vs. immune genes
fisher.test(immunecounts[c("Effector","AllGenes"),c("m7m8select","m7m8notselect")])
#Receptor vs. immune genes
fisher.test(immunecounts[c("Receptor","AllGenes"),c("m7m8select","m7m8notselect")])
#Signaling vs. immune genes
fisher.test(immunecounts[c("Signaling","AllGenes"),c("m7m8select","m7m8notselect")])


#Fisher's Exact test for signficance of different immune classes vs. all immune genes with the m1 vs. m2 model with omega > 2 and with at least 5% of sites
#Effector vs. immune genes
fisher.test(immunecounts[c("Effector","AllGenes"),c("m1m2selectOmega","m1m2notselectOmega")])
#Receptor vs. immune genes
fisher.test(immunecounts[c("Receptor","AllGenes"),c("m1m2selectOmega","m1m2notselectOmega")])
#Signaling vs. immune genes
fisher.test(immunecounts[c("Signaling","AllGenes"),c("m1m2selectOmega","m1m2notselectOmega")])
#None are signficant

#Fisher's Exact test for signficance of different immune classes vs. all immune genes with the m7 vs. m8 model with omega > 2 and with at least 5% of sites
#Effector vs. immune genes
fisher.test(immunecounts[c("Effector","AllGenes"),c("m7m8selectOmega","m7m8notselectOmega")])
#Receptor vs. immune genes
fisher.test(immunecounts[c("Receptor","AllGenes"),c("m7m8selectOmega","m7m8notselectOmega")])
#Signaling vs. immune genes
fisher.test(immunecounts[c("Signaling","AllGenes"),c("m7m8selectOmega","m7m8notselectOmega")])



barplot(immunecounts[c("ImmuneGenes","Receptor","Signaling","Effector"),"prop_m1m2select"],names.arg=c("All Immune(140)","Receptor(50)","Signaling(41)","Effector(46)"), las=2)

immcols <- c("#e7298a","#7570b3","#1b9e77","#d95f02")
names(immcols) <- c("ImmuneGenes","Receptor","Signaling","Effector")


icp <- immunecounts[c("ImmuneGenes","Receptor","Signaling","Effector"),]
icp[,"type"] <- rownames(icp)
icp$type <- factor(icp$type,levels=c("ImmuneGenes","Receptor","Signaling","Effector"))

ggplot(data=icp,aes(type,prop_m1m2select)) + geom_bar(stat="identity",fill=immcols) + geom_hline(yintercept=immunecounts["AllGenes","prop_m1m2select"],linetype=2,size=1.5) + theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1,size=13), axis.title.y = element_text(vjust=1,size=13)) + xlab("Gene Type") + ylab("Proportion Selected") + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) + ylim(0,0.6)

ggplot(data=icp,aes(type,prop_m1m2selectOmega)) + geom_bar(stat="identity",fill= immcols) + geom_hline(yintercept=immunecounts["AllGenes","prop_m1m2selectOmega"],linetype=2,size=1.5) + theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1,size=13), axis.title.y = element_text(vjust=1,size=13)) + xlab("Gene Type") + ylab("Proportion Selected") + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) + ylim(0,0.25)


ggplot(data=icp,aes(type,prop_m7m8select)) + geom_bar(stat="identity",fill= immcols) + geom_hline(yintercept=immunecounts["AllGenes","prop_m7m8select"],linetype=2,size=1.5) + theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1), axis.title.y = element_text(vjust=1)) + xlab("Gene Type") + ylab("Proportion Selected") + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) 


ggplot(data=icp,aes(type,prop_m7m8selectOmega)) + geom_bar(stat="identity",fill= immcols) + geom_hline(yintercept=immunecounts["AllGenes","prop_m7m8selectOmega"],linetype=2,size=1.5) + theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1), axis.title.y = element_text(vjust=1)) + xlab("Gene Type") + ylab("Proportion Selected") + ggtitle("M7 vs. M8") + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())



#Compare omega values across immune classes:
immune_omega <- paml_pval_immune[paml_pval_immune$Class == "Receptor" | paml_pval_immune$Class == "Effector" | paml_pval_immune$Class == "Signaling",c("FDRPval_m1m2","FDRPval_m7m8","Omega_m2","Omega_m8","Class")]
for (i in 1:nrow(immune_omega)){
	if (!is.na(immune_omega[i,"FDRPval_m1m2"]) & immune_omega[i,"FDRPval_m1m2"]>0.05){
		immune_omega[i,"Omega_m2"] <- 1
	}
	if (!is.na(immune_omega[i,"FDRPval_m7m8"]) & immune_omega[i,"FDRPval_m7m8"]>0.05){
		immune_omega[i,"Omega_m8"] <- 1
	}

}

immune_omega_all <- immune_omega
immune_omega_all$Class <- rep("ImmuneGenes",nrow(immune_omega))

immune_omega <- rbind(immune_omega,immune_omega_all)

immune_omega$Class <- factor(immune_omega$Class, levels=c("ImmuneGenes","Receptor","Signaling","Effector"))

ggplot(data=immune_omega, aes(Class,Omega_m2)) + geom_boxplot(fill= immcols) + theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1,size=13), axis.text.y = element_text(size=13),  axis.title = element_text(vjust=1,size=13)) + xlab("Gene Type") + ylab("Omega (dN/dS) Selected") + geom_hline(yintercept=mean(paml_pval_allgenes$Omega_m2,na.rm=T),linetype=2,size=1.5) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) 

ggplot(data=immune_omega, aes(Class,Omega_m8)) + geom_boxplot(fill="darkorchid4") + theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1), axis.title.y = element_text(vjust=1)) + xlab("Gene Type") + ylab("Omega (dN/dS) Selected") + ggtitle("M8 Model") + geom_hline(yintercept=mean(paml_pval_allgenes$Omega_m8,na.rm=T),linetype=2,size=1.5)

for (i in 1:nrow(paml_pval_allgenes)){
	if (!is.na(paml_pval_allgenes[i,"FDRPval_m1m2"]) & paml_pval_allgenes[i,"FDRPval_m1m2"]>0.05){
		paml_pval_allgenes[i,"Omega_m2"] <- 1
	}
	if (!is.na(paml_pval_allgenes[i,"FDRPval_m7m8"]) & paml_pval_allgenes[i,"FDRPval_m7m8"]>0.05){
		paml_pval_allgenes[i,"Omega_m8"] <- 1
	}

}





ks.test(immune_omega[immune_omega$Class=="ImmuneGenes","Omega_m2"],immune_omega[immune_omega$Class=="Receptor","Omega_m2"])
ks.test(immune_omega[immune_omega$Class=="ImmuneGenes","Omega_m2"],immune_omega[immune_omega$Class=="Signaling","Omega_m2"])
ks.test(immune_omega[immune_omega$Class=="ImmuneGenes","Omega_m2"],immune_omega[immune_omega$Class=="Effector","Omega_m2"])

ks.test(immune_omega[immune_omega$Class=="ImmuneGenes","Omega_m8"],immune_omega[immune_omega$Class=="Receptor","Omega_m8"])
ks.test(immune_omega[immune_omega$Class=="ImmuneGenes","Omega_m8"],immune_omega[immune_omega$Class=="Signaling","Omega_m8"])
ks.test(immune_omega[immune_omega$Class=="ImmuneGenes","Omega_m8"],immune_omega[immune_omega$Class=="Effector","Omega_m8"])

#Against all genes
ks.test(paml_pval_allgenes$Omega_m2,immune_omega[immune_omega$Class=="Receptor","Omega_m2"])
ks.test(paml_pval_allgenes$Omega_m2,immune_omega[immune_omega$Class=="Signaling","Omega_m2"])
ks.test(paml_pval_allgenes$Omega_m2,immune_omega[immune_omega$Class=="Effector","Omega_m2"])

ks.test(paml_pval_allgenes$Omega_m8,immune_omega[immune_omega$Class=="Receptor","Omega_m8"])
ks.test(paml_pval_allgenes$Omega_m8,immune_omega[immune_omega$Class=="Signaling","Omega_m8"])
ks.test(paml_pval_allgenes$Omega_m8,immune_omega[immune_omega$Class=="Effector","Omega_m8"])




#Pathway enrichment:

biogenes <- read.table("~/Dropbox/BirdImmuneGeneEvolution/Biosystems/biosystems_gene")
colnames(biogenes) <- c("bsid","GeneID","score")

biopath <- read.table("~/Dropbox/BirdImmuneGeneEvolution/Biosystems/bsid2info",sep="\t", fill=T,header=F)

#Assign paml results to pathways from biosystems.
paml_pval_innatedb_path <- merge(biogenes,paml_pval_innatedb,by.x="GeneID",by.y="entrezgene")
paml_pval_allgenes_path <- merge(biogenes,paml_pval_allgenes_ncbi,by.x="GeneID",by.y="entrezgene")

#How many of our paml results could be assigned ot pathways?
length(unique(paml_pval_allgenes_path$GeneID))
#3581
#Of how many?
length(unique(paml_pval_allgenes_ncbi$entrezgene))
#11574
#setdiff(unique(paml_pval_allgenes_ncbi$entrezgene),unique(paml_pval_allgenes_path$GeneID))

#Seems to be a number of genes missing?




#Try for KEGG pathway enrichment
require(DOSE)
require(clusterProfiler)

m2_signames <- paml_pval_allgenes_ncbi[paml_pval_allgenes_ncbi$FDRPval_m1m2<0.05,"entrezgene"]

m2_k <- enrichKEGG(m2_signames,organism="gga",pvalueCutoff=0.05,pAdjustMethod="BH",qvalueCutoff=0.01)

m2_sigOmeganames <- paml_pval_allgenes_ncbi[paml_pval_allgenes_ncbi$FDRPval_m1m2<0.05 & paml_pval_allgenes_ncbi$Omega_m2>2 & paml_pval_allgenes_ncbi$Prop_m2>0.05,"entrezgene"]

m2_omega_k <- enrichKEGG(m2_sigOmeganames,organism="gga",pvalueCutoff=0.05,pAdjustMethod="BH",qvalueCutoff=0.05)

summary(m2_omega_k)



#Visualize pathways
require(pathview)

#Extract omegas from m2 model and assign non-sig to 0.01 (for visualization purposes)
rownames(paml_pval_allgenes_ncbi) <- paml_pval_allgenes_ncbi$entrezgene
omega_m2 <- paml_pval_allgenes_ncbi[,"Omega_m2"]
omega_m2[paml_pval_allgenes_ncbi$FDRPval_m1m2>0.05] <- 1
names(omega_m2) <- rownames(paml_pval_allgenes_ncbi)

pv.out.04020 <- pathview(gene.data=omega_m2,pathway.id="04620",species="gga",both.dirs=list(gene=F),mid="#e5f5f9",high="#2ca25f",limit=list(gene=5),kegg.native=f,key.pos="topleft",out.suffix="TLR_Signaling")

pv.out.05164 <- pathview(gene.data=omega_m2,pathway.id="05164",species="gga",both.dirs=list(gene=F),mid="#e5f5f9",high="#2ca25f",limit=list(gene=5),kegg.native=F,key.pos="topright",out.suffix="Influenza_A")

pv.out.04514 <- pathview(gene.data=omega_m2,pathway.id="04514",species="gga",both.dirs=list(gene=F),mid="#e5f5f9",high="#2ca25f",limit=list(gene=5),kegg.native=F,key.pos="topright",out.suffix="CAMs")

pv.out.05168 <- pathview(gene.data=omega_m2,pathway.id="05168",species="gga",both.dirs=list(gene=F),mid="#e5f5f9",high="#2ca25f",limit=list(gene=5),kegg.native=T,key.pos="topright",out.suffix="HerpesSimplexInfection")

pv.out.04060 <- pathview(gene.data=omega_m2,pathway.id="04060",species="gga",both.dirs=list(gene=F),mid="#e5f5f9",high="#2ca25f",limit=list(gene=5),kegg.native=T,key.pos="topright",out.suffix="CytokineCytokineReceptorInteraction")


#pv.out <- pathview(gene.data=omega_m2,pathway.id="04620",species="gga",both.dirs=list(gene=F),mid="#7fbf7b",high="#af8dc3",limit=list(gene=5),na.col="gray",kegg.native=F)

omega_m2_sig <- paml_pval_allgenes_ncbi[,c("FDRPval_m1m2","Omega_m2")]
m2_sig <- rep(0,nrow(paml_pval_allgenes_ncbi))
m2_sig[paml_pval_allgenes_ncbi$FDRPval_m1m2<0.05] <- 1
names(m2_sig) <- rownames(paml_pval_allgenes_ncbi)
m2_sig <- names(m2_sig[m2_sig==1])


paml_pval_allgenes_ncbi[paml_pval_allgenes_ncbi$Omega_m2>1 & paml_pval_allgenes_ncbi$FDRPval_m1m2>0.05,]


pv.out <- pathview(gene.data=omega_m2_sig[,1],pathway.id="04620",species="gga",both.dirs=list(gene=F),mid="#e5f5f9",high="#2ca25f",limit=list(gene=15),na.col="gray",kegg.native=F,out.suffix="Omega_m2_sig")

pv.out <- pathview(gene.data=m2_sig,pathway.id="04620",species="gga",both.dirs=list(gene=F),mid="#e5f5f9",high="#2ca25f",na.col="gray",kegg.native=F,out.suffix="m2_sig",discrete=list(gene=T))



require(ape)

tree <- read.tree("~/Dropbox/BirdImmuneGeneEvolution/11700.final_spt.nwk")
plot(tree,label.offset=0.7,cex=0.8,no.margin=T)

sp <- read.csv("~/Dropbox/BirdImmuneGeneEvolution/species_list.csv",header=T)
rownames(sp) <- sp$Short.name

old_label <- tree$tip.label
tree$tip.label <- as.character(sp[tree$tip.label,"Binomial.name"])
plot(tree,label.offset=0.7,cex=0.8,no.margin=T)

tree <- read.tree("~/Dropbox/BirdImmuneGeneEvolution/11700.final_spt.nwk")
tree$tip.label <- as.character(sp[tree$tip.label,"Common.name"])
plot(tree,label.offset=0.7,cex=0.8,no.margin=T)


