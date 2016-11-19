###########Which HOGs don't have PAML results?
paml_res <- read.table("~/Dropbox/BirdImmuneGeneEvolution/PAML/paml_res_allhogs_rewrite.txt",header=T)
hog_ids_res <- unique(paml_res$hog)

#Original list
all_hogs <- read.table("~/Dropbox/BirdImmuneGeneEvolution/PAML/all_hogs")
all_hogs <- all_hogs$V1

#See how many results each model has
hog_ids_res_nums <- table(paml_res$hog)

#Collect all those that do not have 4
hog_rerun <- names(hog_ids_res_nums[hog_ids_res_nums!=4])

#Which hogs don't have any results?
rerun2 <- setdiff(all_hogs,hog_ids_res)

#Combine the lists
all_rerun <- c(hog_rerun,as.character(rerun2))
length(all_rerun)
#1029 HOGs need to be rerun.


all_rerun <- data.frame(all_rerun)

write.table(all_rerun,"~/Dropbox/BirdImmuneGeneEvolution/PAML/hogs_to_rerun.txt",sep="\t",row.names=F,col.names=F,quote=F)

#Which of the 1029 need to be rerun?
all_rerun <- read.table("~/Dropbox/BirdImmuneGeneEvolution/PAML/hogs_to_rerun.txt")
all_rerun <- all_rerun$V1

all_rerun_res <- read.table("~/Dropbox/BirdImmuneGeneEvolution/PAML/hogs_to_rerun_results.txt",header=T)
all_rerun_res_m7 <- read.table("~/Dropbox/BirdImmuneGeneEvolution/PAML/hogs_to_rerun_results_m7Only.txt",header=T)
rownames(all_rerun_res_m7) <- all_rerun_res_m7$hog
all_rerun_res_m8 <- read.table("~/Dropbox/BirdImmuneGeneEvolution/PAML/hogs_to_rerun_results_m8Only.txt",header=T)
rownames(all_rerun_res_m8) <- all_rerun_res_m8$hog

all_rerun_res$omega <- as.character(all_rerun_res$omega)
all_rerun_res_m7$omega <- as.character(all_rerun_res_m7$omega)
all_rerun_res_m8$omega <- as.character(all_rerun_res_m8$omega)

need_m1 <- vector()
need_m2 <- vector()
need_m7 <- vector()
need_m8 <- vector()

#Integrate separate m7 and m8 results into full results table - only add if missing.
for (i in 1:length(all_rerun)){
	hog <- as.character(all_rerun[i])
	
	if (!("1" %in% all_rerun_res[all_rerun_res$hog==hog,"model_num"])){
		need_m1[(length(need_m1)+1)] <- hog}
	else {}
	
	if (!("2" %in% all_rerun_res[all_rerun_res$hog==hog,"model_num"])){
		need_m2[(length(need_m2)+1)] <- hog}
	else {}
	
	if (!("7" %in% all_rerun_res[all_rerun_res$hog==hog,"model_num"])){
		if (hog %in% all_rerun_res_m7$hog){
			all_rerun_res[(nrow(all_rerun_res)+1),] <- all_rerun_res_m7[as.character(hog),]
			print(hog)
		}
		else {need_m7[(length(need_m7)+1)] <- hog}
		}
	else {}

	if (!("8" %in% all_rerun_res[all_rerun_res$hog==hog,"model_num"])){
		if (hog %in% all_rerun_res_m8$hog){
			all_rerun_res[(nrow(all_rerun_res)+1),] <- all_rerun_res_m8[as.character(hog),]
			print(hog)
		}
		else {need_m8[(length(need_m8)+1)] <- hog}
		}
	else {}

}

write.table(need_m1,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/hogs_to_rerun_m1",sep="\t",row.names=F,col.names=F,quote=F)
write.table(need_m2,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/hogs_to_rerun_m2",sep="\t",row.names=F,col.names=F,quote=F)
write.table(need_m7,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/hogs_to_rerun_m7",sep="\t",row.names=F,col.names=F,quote=F)
write.table(need_m8,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/hogs_to_rerun_m8",sep="\t",row.names=F,col.names=F,quote=F)



#Checking which to rerun with models 2a and 8a.
all_hogs <- read.table("~/Dropbox/BirdImmuneGeneEvolution/PAML/all_hogs")
all_hogs <- as.character(all_hogs$V1)

all_res_m2a <- read.table("~/Dropbox/BirdImmuneGeneEvolution/PAML/allhogs_results_m2aOnly.txt",header=T)
rownames(all_res_m2a) <- as.character(all_res_m2a$hog)

all_res_m8a <- read.table("~/Dropbox/BirdImmuneGeneEvolution/PAML/allhogs_results_m8aOnly.txt",header=T)
rownames(all_res_m8a) <- all_res_m8a$hog

need_m2a <- setdiff(all_hogs,rownames(all_res_m2a))
need_m8a <- setdiff(all_hogs,rownames(all_res_m8a))

write.table(need_m2a,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/hogs_to_rerun_m2a",sep="\t",row.names=F,col.names=F,quote=F)
write.table(need_m8a,file="~/Dropbox/BirdImmuneGeneEvolution/PAML/hogs_to_rerun_m8a",sep="\t",row.names=F,col.names=F,quote=F)


