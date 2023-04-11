#library(rhmmer)
library(dplyr)
library(stringr)
source('src/rhmmer_functions.R')
#arg[1] is the dbcan annotation
#arg[2] is the pfam annotation
#arg[3] is the dbcan family info
#arg[4] is the pfam family info
#arg[5] is the output table with both annotations


args = commandArgs(trailingOnly=TRUE)


###Extract top hits function
top_hits <- function(table){
  genes <- unique(table$domain_name)
  top_hits <- table[0,]
  for (i in 1:length(genes)){
    subtable <- table[table$domain_name == genes[i],]
    min_evalue <- subtable[which.min(subtable$sequence_evalue),]
    top_hits <- rbind(top_hits,min_evalue)
  }
  
  return(top_hits)
}


##DBCAN annotations
dbcan <- read_tblout(args[1])
dbcan <- top_hits(dbcan)
dbcan$query_name <- sapply(strsplit(dbcan$query_name,"_"), `[`, 1) 
dbcan$query_name <- str_remove(dbcan$query_name, '.hmm')

dbcan_description <- read.delim(args[3], header = F)
colnames(dbcan_description)[1]<- 'dbcan_family'
colnames(dbcan_description)[2]<- 'dbcan_description'


dbcan_final <- merge(dbcan[,c(1,3)], dbcan_description, by.x = 'query_name', by.y = 'dbcan_family', all.x=T)
colnames(dbcan_final)[1]<- 'dbcanid'
colnames(dbcan_final)[2]<- 'gene'
##Group different annotations of same gene
dbcan_final <- dbcan_final %>% group_by_at(vars(gene)) %>%
  summarize_all(paste, collapse=",")
###Check problems of family names with a _

##PFAM annotations
pfam <- read_tblout(args[2])
pfam <- top_hits(pfam)
pfam <- pfam[,c(1,4)]
pfam$query_accession <- sapply(strsplit(pfam$query_accession,"[.]"), `[`, 1)

pfam_description <- read.delim2(args[4], header = F)
pfam_description <-pfam_description[,c(1,5)]
colnames(pfam_description)[1] <- 'query_accession'
colnames(pfam_description)[2] <- 'pfam_description'

pfam_final <- merge(pfam, pfam_description, by = 'query_accession', all.x = T)
colnames(pfam_final)[1]<- 'pfamid'
colnames(pfam_final)[2]<- 'gene'
##Group different annotations of same gene
pfam_final <- pfam_final %>% group_by_at(vars(gene)) %>%
  summarize_all(paste, collapse=",")


###Merge both annotations in one table
hmm_annotations <- merge(pfam_final, dbcan_final, by= 'gene', all = T)
write.table(hmm_annotations, args[5], row.names=FALSE)