####################################
###     Guilermo Guerrero        ###
### G.GuerreroEgido@nioo.knaw.nl ###
###           2021               ###
###     Leiden University        ###
###                              ###
####################################

##Install packages if not installed
listOfPackages <- c("shiny","readr","UpSetR","data.table","dplyr",
                    "dbscan","shinydashboard","scales","shinythemes",
                    "plotly","DT", "reshape2","ggplot2","gridExtra","gplots",
                   "gtools","RColorBrewer","viridis","rmarkdown", "ape",
                   "heatmaply","pals","stringdist","tidyr",
                   "stringr","parallelDist","hexbin", "parallel", "BiocManager")
for (i in listOfPackages){
     if(! i %in% installed.packages()){
         install.packages(i, dependencies = TRUE)
     }
}

##Install bioconductor packages if not installed

listOfBCPackages <- c("Biostrings","BiocGenerics")

for (i in listOfBCPackages){
     if(! i %in% installed.packages()){
         BiocManager::install(i)
     }
}


library(shiny)
library(readr)
library(UpSetR)
library(data.table)
library(dplyr)
library(dbscan)
library(shinydashboard)
library(scales)
library(shinythemes)
library(plotly)
library(DT)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(gplots)
library(gtools)
library(RColorBrewer)
library(viridis)
library(rmarkdown)
library(heatmaply)
library(pals)
library(ape) #Bioconductor manager
library(stringdist)
library(BiocGenerics) #Bioconductor manager
library(Biostrings) #Bioconductor manager
library(tidyr)
library(stringr)
library(parallel)
library(parallelDist)
library(hexbin)

###COREPAN ANALYSIS###
###                ###

####FULL MATRIX UPSET plots
source('src/corepan_functions.R')
source('src/stats_functions.R')
source('src/text_shiny.R')


###Inputs
mapping_file <-  read.table('input/mapping_file.txt', header = T)
matrix <- read_delim('input/MEGAMATRIX_renamed.txt', col_names = T, quote = "\"")
all_proteins <- ape::read.FASTA('input/combined_proteins.fasta', type="AA")
big_scape_matrix <- read.table('input/big_scape_binary_table_renamed.txt', header = T, check.names = F)
gcf_annotation <- read.csv('input/annotation.txt', header = T)
names_equivalence = read.table('input/names_equivalence.txt', header = T)


###
mapping_file$Species <- sapply(strsplit(mapping_file$Sample,"_"), `[`, 2)
matrix = as.data.frame(matrix)
n_samples <- grep("completeness", colnames(matrix)) - 1
matrix[,2:n_samples] <- sapply(matrix[,2:n_samples],as.numeric)



###KEGG MATRIX UPSET PLOTS
kegg_matrix <- matrix[,c(2:n_samples, (n_samples + 4))]
kegg_agg_matrix = aggregate(kegg_matrix[,1:(n_samples - 1)], by = list(id = kegg_matrix[,n_samples]), FUN=sum)

###COG MATRIX UPSET PLOTS
cog_matrix <- matrix[,c(2:n_samples, (n_samples + 6))]
cog_matrix <- data.table(cog_matrix)
cog_matrix[is.na(cog_matrix$cogid)]$cogid <- 'Unknown'
cog_matrix <- as.data.frame(cog_matrix)
cog_agg_matrix <- aggregate(cog_matrix[,1:(n_samples - 1)], by = list(id = cog_matrix[,n_samples]), FUN=sum)


###PFAM MATRIX UPSET PLOTS
pfam_matrix <- matrix[,c(2:n_samples, (n_samples + 8))]
pfam_agg_matrix <- aggregate(pfam_matrix[,1:(n_samples - 1)], by = list(id = pfam_matrix[,n_samples]), FUN=sum)


###PROKKA MATRIX UPSET PLOTS
prokka_matrix <- matrix[,c(2:n_samples, (n_samples + 12))]
prokka_agg_matrix <- aggregate(prokka_matrix[,1:(n_samples - 1)], by = list(id = prokka_matrix[,n_samples]), FUN=sum)



###DBCAN MATRIX UPSET PLOTS
dbcan_matrix <- matrix[,c(2:n_samples, (n_samples + 10))]
dbcan_agg_matrix <- aggregate(dbcan_matrix[,1:(n_samples - 1)], by = list(id = dbcan_matrix[,n_samples]), FUN=sum)



#Load cog annotations and clean it
cog_extended_annotations <- read.csv('src/cog_annotation_groups.csv', header = F)
cog_extended_annotations <- cog_extended_annotations[,c(1,2,8)]
colnames(cog_extended_annotations)[1] <- 'id'
colnames(cog_extended_annotations)[3] <- 'general processes'
colnames(cog_extended_annotations)[2] <- 'processes'

cog_processes <- cog_options(cog_agg_matrix, cog_extended_annotations)
df_colors <- data.frame(processes = sort(unique(cog_extended_annotations$processes)), color = assign_colors(sort(unique(cog_extended_annotations$processes))))


####BGCs 
GCF_matrix <- process_bgc_table(big_scape_matrix,gcf_annotation)
mapping_file_GCF <- mapping_file[mapping_file$Lifestyle != 'New',]




dir.create('corepan_analysis')
####Extract core_genes
core_genes_out <- core_genes_extraction(matrix, n_samples)
core_genes <- core_genes_out[[1]]
write.table(core_genes,'corepan_analysis/core_genes.txt', row.names = F)
non_core_genes <- core_genes_out[[2]]

##COG hierarchy of coregenes
cog_core_genes <- all_core_genes_cog(core_genes, cog_extended_annotations)




##Extract singletons genes
extract_singleton_genes <- extract_singletons(non_core_genes)
singletons <- extract_singleton_genes[[1]]
no_singletons <- extract_singleton_genes[[2]]

##Remove genes with more than 90 % of 0s
clean_matrix <- remove_genes_with_0s(no_singletons)
write.table(clean_matrix,'corepan_analysis/clean_matrix.txt', row.names = F)

##Extract singletons for each sample(output automatic in directory singletons_per_sample/ )
singletons_table <- singleton_per_sample(singletons)


###Exploratory ANALYSIS###
###                    ###

##Full matrix



n_samples <- grep("completeness", colnames(clean_matrix)) - 1
clean_matrix[2:n_samples] <- apply(clean_matrix[2:n_samples], 2, function(x) ifelse(x > 1, 1, x))

full_matrix <- clean_matrix[,1:n_samples] #include number depending on number of samples
#remove rows where all columns are the same
full_matrix<- remove_constant_columns(full_matrix)

rownames(full_matrix) <- full_matrix$clusters 
full_matrix$clusters <- NULL
full_matrix <- as.matrix(full_matrix)

full_pca_input <- pcoa_function(full_matrix)


#full_pca_input <- pcoa_function(full_matrix)
full_pca <- merge(full_pca_input[[1]], mapping_file, by = 'Sample', all.x = T)
full_eig <- full_pca_input[[2]]

full_dend <- dendogram(full_matrix)



###GCF matrix
gcf_pca_input <- GCF_matrix
gcf_pca_input[2:(ncol(gcf_pca_input)-1)] <- apply(gcf_pca_input[2:(ncol(gcf_pca_input)-1)], 2, function(x) ifelse(x > 1, 1, x))
rownames(gcf_pca_input) <-gcf_pca_input$id
gcf_pca_input$id <- NULL
gcf_pca_input$BGC.Class <- NULL
gcf_pca_input <- as.matrix(gcf_pca_input)
GCF_pca_input <- pcoa_function(gcf_pca_input)
GCF_pca <- merge(GCF_pca_input[[1]], mapping_file, by = 'Sample', all.x = T)
GCF_eig <- GCF_pca_input[[2]]

GCF_dend <- dendogram(gcf_pca_input)


###KEGG matrix ###
kegg_matrix <- matrix[,c(2:n_samples, grep("keggid", colnames(matrix)))]
kegg_matrix <- aggregate_matrix(kegg_matrix, 'kegg', mapping_file, n_samples)
kegg_matrix[kegg_matrix > 1] <- 1
kegg_pca_input <- pcoa_function(kegg_matrix)
kegg_pca <- merge(kegg_pca_input[1], mapping_file, by = 'Sample', all.x = T)
kegg_eig <- kegg_pca_input[[2]]


kegg_dend <- dendogram(kegg_matrix)




###COG matrix ###
cog_matrix <- matrix[,c(2:n_samples, grep("cogid", colnames(matrix)))]
cog_matrix <- aggregate_matrix(cog_matrix, 'cog', mapping_file, n_samples)
cog_matrix[cog_matrix > 1] <- 1
cog_pca_input <- pcoa_function(cog_matrix) 
cog_pca <- merge(cog_pca_input[[1]], mapping_file, by = 'Sample', all.x = T)
cog_eig <- cog_pca_input[[2]]

cog_dend <- dendogram(cog_matrix)

###PFAM matrix ###
pfam_matrix <- matrix[,c(2:n_samples, grep("pfamid", colnames(matrix)))]
pfam_matrix <- aggregate_matrix(pfam_matrix, 'pfam', mapping_file, n_samples)
pfam_matrix[pfam_matrix > 1] <- 1
pfam_pca_input <- pcoa_function(pfam_matrix) 
pfam_pca <- merge(pfam_pca_input[[1]], mapping_file, by = 'Sample', all.x = T)
pfam_eig <- pfam_pca_input[[2]]

pfam_dend <- dendogram(pfam_matrix)

###PROKKA matrix ###
prokka_matrix <- matrix[,c(2:n_samples, grep("prokkadescription", colnames(matrix)))]
prokka_matrix <- prokka_matrix[prokka_matrix$prokkadescription != 'hypothetical protein',] 
prokka_matrix <- aggregate_matrix(prokka_matrix, 'prokka', mapping_file, n_samples)
prokka_matrix[prokka_matrix > 1] <- 1
prokka_pca_input <- pcoa_function(prokka_matrix) 
prokka_pca <- merge(prokka_pca_input[[1]], mapping_file, by = 'Sample', all.x = T)
prokka_eig <- prokka_pca_input[[2]]

prokka_dend <- dendogram(prokka_matrix)

###DBCAN matrix ###
dbcan_matrix <- matrix[,c(2:n_samples, grep("dbcanid", colnames(matrix)))]
dbcan_matrix <- aggregate_matrix(dbcan_matrix, 'dbcan', mapping_file, n_samples)
dbcan_matrix[dbcan_matrix > 1] <- 1
dbcan_pca_input <- pcoa_function(dbcan_matrix)
dbcan_pca <- merge(dbcan_pca_input[[1]], mapping_file, by = 'Sample', all.x = T)
dbcan_eig <- dbcan_pca_input[[2]]

dbcan_dend <- dendogram(dbcan_matrix)




mapping_options <- as.list(c(colnames(mapping_file[,2:ncol(mapping_file)]),'clusters_knn', 'clusters_hdbscan', 'clusters_dendogram' ))
mapping_options_PCA <- as.list(colnames(mapping_file[,2:ncol(mapping_file)]))

cog_process_options <- as.list(c('All', cog_processes ))



###Free RAM
rm()

