library(data.table)
library(dplyr)


args = commandArgs(trailingOnly=TRUE)
##arg[1] is the intermediate_files/matrices_cdhit/cd-hit/orthogous_groups.txt
##arg[2] is the 'MEGAMATRIX.txt'
##arg[3] is the cogg annotatuion intermediate_files/cog_annotation/COG_annotation_clean.faa.finalcog
##arg[4] is the kegg annotation 'intermediate_files/kegg_annotation/KEGG_annotation_clean.faa.finalkegg'
##arg[5] is the output name of matrix with all annotations







###
###FIX gene and description column of new samples in MEGAMATRIX
###
orth_file <- read.delim2(args[1], header = F)
orth_file$V3 <- sapply(strsplit(orth_file$V3,";"), `[`, 1)
orth_file$cluster_n <- as.numeric( sapply(strsplit(orth_file$V1,"_"), `[`, 2))
orth_file <- orth_file[order(orth_file$cluster_n),]


MEGAMATRIX <- read.table(args[2], header = T)
MEGAMATRIX$cluster_n <- as.numeric(sapply(strsplit(MEGAMATRIX$clusters,"_"), `[`, 2))
MEGAMATRIX <- MEGAMATRIX[order(MEGAMATRIX$cluster_n),]


MEGAMATRIX$descriptions <- orth_file$V2
MEGAMATRIX <- data.table(MEGAMATRIX)
orth_file <- data.table(orth_file)


MEGAMATRIX[is.na(MEGAMATRIX$gene),]$gene <- sapply(strsplit(MEGAMATRIX[is.na(MEGAMATRIX$gene),]$descriptions,";"), `[`, 1)


###Add the annotation of new clusters COG + KEGG

cog <- read.table(args[3], header = T)
colnames(cog)[2:3] <- paste0(colnames(cog)[2:3], '.y')
cog <- data.table(cog)


mi_fun <- function(i, a) {
   new_cog <- a[a %in% strsplit(i,";")[[1]]]
   if (length(new_cog) > 0){
     return_value <- new_cog
   } 
   else{
     return_value <- NA
   }
   
  return(return_value)
}


MEGAMATRIX$Add_cog <- sapply(MEGAMATRIX$descriptions, mi_fun, a = cog$gene)
MEGAMATRIX <- merge(MEGAMATRIX,   cog, by.x = 'Add_cog', by.y = 'gene', all.x = T)
MEGAMATRIX[!is.na(MEGAMATRIX$Add_cog)]$cogid<- MEGAMATRIX[!is.na(MEGAMATRIX$Add_cog)]$cogid.y
MEGAMATRIX[!is.na(MEGAMATRIX$Add_cog)]$cog_description <- MEGAMATRIX[!is.na(MEGAMATRIX$Add_cog)]$cog_description.y

MEGAMATRIX$cogid.y <- NULL
MEGAMATRIX$cog_description.y <- NULL
MEGAMATRIX$Add_cog <- NULL



###KEGG
kegg <- read.table(args[4], header = T)
colnames(kegg)[2:3] <- paste0(colnames(kegg)[2:3], '.y')
kegg <- data.table(kegg)


MEGAMATRIX$Add_kegg <- sapply(MEGAMATRIX$descriptions, mi_fun, a = kegg$gene)
MEGAMATRIX <- merge(MEGAMATRIX,   kegg, by.x = 'Add_kegg', by.y = 'gene', all.x = T)
MEGAMATRIX[!is.na(MEGAMATRIX$Add_kegg)]$keggid<- MEGAMATRIX[!is.na(MEGAMATRIX$Add_kegg)]$keggid.y
MEGAMATRIX[!is.na(MEGAMATRIX$Add_kegg)]$kegg_description <- MEGAMATRIX[!is.na(MEGAMATRIX$Add_kegg)]$kegg_description.y

MEGAMATRIX$keggid.y <- NULL
MEGAMATRIX$kegg_description.y <- NULL
MEGAMATRIX$Add_kegg <- NULL



write.table(MEGAMATRIX, args[5], row.names = F)

