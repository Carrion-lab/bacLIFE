library(igraph)
library(stringr)

###remove redundant genomes:

mash <- read.table('mash/combined_results_filtered', sep = '\t', header = F)
mash <-  mash[mash$V1!=mash$V2,]
mash <- mash[mash$V3 < 0.01,]

ig <- graph_from_data_frame(mash[,1:3], directed = F)
cl <- components(ig)

clusters <- as.data.frame(cl$membership)
clusters$Sample <- rownames(clusters)

clusters$Sample <- str_remove(clusters$Sample, 'genomes/')


lista_all_files <- data.frame( Sample =  list.files('genomes/', pattern = '.fa'), downloaded = 'downloaded' )


M <- merge(clusters, lista_all_files, by = 'Sample', all = T)
uniques <- M[is.na(M$`cl$membership`),]
uniques <- uniques$Sample

clusters_chosen = clusters[!duplicated(clusters$`cl$membership`),]

to_use <- c(uniques, clusters_chosen$Sample)


###Generate metadata table of mash clusters

M <- data.table::data.table(M)
M$`cl$membership` <- as.numeric(M$`cl$membership`)
M <- M[order(M$`cl$membership`),]
n_clusters <- max(M[!is.na(M$`cl$membership`)]$`cl$membership`) +1 
n_na <- nrow(M[is.na(M$`cl$membership`)])
M[is.na(M$`cl$membership`)]$`cl$membership` <- seq(n_clusters, (n_clusters + n_na - 1))

M$downloaded <- NULL
M$representative <- 'clustered'
M[M$Sample %in% to_use]$representative <- 'representative'
colnames(M)[2] <- 'Cluster_membership'
write.table(M, 'cluster_membership.txt', row.names = F)



##Move files to rename

file.copy(paste0('genomes/', to_use), 'genomes_renamed/')
'Files were renamed and copied in genomes_renamed/'




