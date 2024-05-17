library(reshape2)
library(tidyverse)
library(stringr)
library(dplyr)
library(readr)
library(Biostrings)
library(data.table)
library(parallel)

cl <- makeCluster(10)
args = commandArgs(trailingOnly=TRUE)

tsv = data.table(read_table(args[1], col_names = F))
colnames(tsv) <- c('V1', 'V2')

n_clusters_mmseq = length(unique(tsv$V1))


tsv2 <- tsv[, .( V2 = str_c(V2, collapse = ',') ), by = V1]

#write.table(tsv2,'diamond/mmseq2_clusters.txt', row.names = F, quote = F)


'Open MCL output, collapse it and clean it\n'

mcl <- data.table(read_csv( args[2], col_names = F))

#Load unaligned headers 
unalin <- read.delim2(args[3], header = F)
unalin <- data.frame(X1 = str_remove_all(unalin$V1, '>'))

mcl <- dplyr::bind_rows(mcl, unalin)


#Name clusters

mcl$clusters <- seq(nrow(mcl))
mcl$clusters <- sprintf("cluster_%06d", mcl$clusters)


mcl_melt <- mcl[, c(X1 = strsplit(X1, "\t")), by = clusters]

'Join mmseq2 clusters with MCL clusters\n'

M <- merge(mcl_melt, tsv2, by.x= 'X1', by.y = 'V1')
M <- M[,c(2,3)]
colnames(M)[2] <- 'descriptions'
M <- data.table(M)
M <- M[, .( descriptions = str_c(descriptions, collapse = ',') ), by = clusters]
M <- M[order(clusters)]
M$gene <- sapply(strsplit(M$descriptions,","), `[`, 1)

#write.table(M, 'clusters_info.txt', row.names = F, sep = '\t', quote = F)

'Create 1/0 matrix\n' 

bacteria_name <- function(string){
  b <- stringr::str_split(string, ',')[[1]]
  c <- sapply(stringr::str_split(b,"[|]"), `[`, 1)
  d <- stringr::str_c(c, collapse = ",")
  return(d)
}

MEGAMATRIX <- M
MEGAMATRIX$descriptions <- parSapply(cl, MEGAMATRIX$descriptions , FUN = bacteria_name)
MEGAMATRIX$gene <- NULL



tbl <- table(splitstackshape:::cSplit(MEGAMATRIX, "descriptions", sep = ",", direction = "long"))
tbl <- data.table(as.data.frame(tbl))
tbl <-dcast.data.table(tbl, clusters ~ descriptions, value.var = 'Freq')

#lista <- as.list(MEGAMATRIX$descriptions)
#binary_matrix <- as.data.frame((splitstackshape:::charMat(listOfValues = lista, fill = 0L)))
#binary_matrix <- cbind(clusters = MEGAMATRIX$clusters, binary_matrix)
binary_matrix <- merge(tbl, M, by = 'clusters')



####Extract sequence of representative proteins
#system("bioawk -c fastx '{ print $name, length($seq) }' < intermediate_files/combined_proteins/combined_proteins.fasta >intermediate_files/combined_proteins/length_genes.txt")
gene_lengths <- read.table('intermediate_files/combined_proteins/length_genes.txt', header = F)

##Function to extract longest gene

extract_long_genes <- function(character, gene_length){
  a <- as.character(stringr::str_split(character, ',')[[1]])
  sub <- data.table::data.table(gene_length[gene_length$V1 %in% a,])
  longest_gene <- sub[which.max(sub$V2)]$V1
  return(longest_gene)
}

'Obtain representative genes\n'

genes_description <- binary_matrix$descriptions

new_rep <- parSapply(cl, genes_description, extract_long_genes, gene_length = gene_lengths)

binary_matrix$gene <- new_rep

write.table(binary_matrix, args[4], row.names = F)

'Obtain representative sequences\n'

representative_genes <- binary_matrix$gene
fasta <-  Biostrings::readAAStringSet(args[5])
subset <- fasta[representative_genes]
writeXStringSet(subset, args[6])




