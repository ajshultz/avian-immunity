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