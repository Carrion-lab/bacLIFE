##############################################
###         Guilermo Guerrero              ###
### g.guerrero.egido@biology.leidenuniv.nl ###
###               2021                     ###
###          Leiden University             ###
###                                        ###
##############################################


##This script takes as input the MEGAMATRIX and the mapping file to test its predictability using random forest
##User can provide its own metadata as a tab separated file or use one of the metadata files we provide (see examples to see format) 
##This script outputs the ROC plots for the predictions of each class in the metadata 
##If good model performance is observed in the ROC plots, metadata can be augmented for genomes labeled as 'Unknown' with the script 'classifier_augment.R'



library(stringr)
library(dplyr)
library(readr)
library(data.table)

####Functions

`%!in%` = Negate(`%in%`)

core_genes_extraction <- function(matrix, n_samples){
  zeros <- rowSums(matrix[2:n_samples] == 0)
  matrix$n_zeros <- zeros
  core_genes <- matrix[matrix$n_zeros < (n_samples - (n_samples*0.9)),]
  non_core_genes <- matrix[matrix$n_zeros > (n_samples - (n_samples*0.9)),]
  core_genes <- core_genes[1:(length(core_genes)-1)]
  return_list <- list(core_genes, non_core_genes)
  return(return_list)
}


extract_singletons <- function(matrix){
  n_samples <- grep("completeness", colnames(matrix)) - 1
  zeros <- rowSums(matrix[2:n_samples] == 0)
  matrix$n_zeros <- zeros
  singletons <- matrix[matrix$n_zeros == (n_samples - 2),]
  singletons <- singletons[1:(length(singletons)-1)]
  no_singletons <- matrix[!(matrix$n_zeros == (n_samples - 2)),]
  no_singletons <- no_singletons[1:(length(no_singletons)-1)]
  
  return_list <- list(singletons, no_singletons)
  return(return_list)
  
}


####Function extract genes with more than 90 % of 0s
remove_genes_with_0s <- function(matrix){
  n_samples <- grep("completeness", colnames(matrix)) - 1
  zeros <- rowSums(matrix[2:n_samples] == 0)
  matrix$n_zeros <- zeros
  no_singletons <- matrix[matrix$n_zeros < (n_samples - (n_samples*0.01)),]
  no_singletons <- no_singletons[1:(length(no_singletons)-1)]
  
  return(no_singletons)
  
}


remove_constant_columns <- function(matrix){
  n_samples <- ncol(matrix)
  keep <- apply(matrix[2:n_samples], 1, function(x) length(unique(x[!is.na(x)])) != 1)
  new_matrix <-  matrix[keep,]
  return(new_matrix)
}





args = commandArgs(trailingOnly=TRUE)

###Load metadata
mapping_file <- read.table(args[1], header = T)
column2use <- args[3]
mapping_file <- mapping_file[,c('Sample', column2use)]
colnames(mapping_file)[2] <- 'Lifestyle'
###load and filter matrix
matrix <- read_delim(args[2], col_names = T, quote = "\"")
matrix = as.data.frame(matrix)
n_samples <- grep("completeness", colnames(matrix)) - 1
matrix[,2:n_samples] <- sapply(matrix[,2:n_samples],as.numeric)

#Extract core_genes
dir.create('corepan_analysis', showWarnings = F)
core_genes_out <- core_genes_extraction(matrix, n_samples)
core_genes <- core_genes_out[[1]]
non_core_genes <- core_genes_out[[2]]

#Extract singletons genes
extract_singleton_genes <- extract_singletons(non_core_genes)
singletons <- extract_singleton_genes[[1]]
no_singletons <- extract_singleton_genes[[2]]

#Remove genes with more than 95% of 0s
clean_matrix <- remove_genes_with_0s(no_singletons)


#Prepare filtered matrix
n_samples <- grep("completeness", colnames(clean_matrix)) - 1
clean_matrix[2:n_samples] <- apply(clean_matrix[2:n_samples], 2, function(x) ifelse(x > 1, 1, x))
full_matrix <- clean_matrix[,1:n_samples]
full_matrix<- remove_constant_columns(full_matrix)

rownames(full_matrix) <- full_matrix$clusters 
full_matrix$clusters <- NULL
full_matrix <- as.matrix(full_matrix)




mtrx = t(full_matrix)
mtrx = as.data.frame(mtrx)
mtrx$Sample = rownames(mtrx)
MMM = merge(mtrx, mapping_file, by = 'Sample')
MMM = data.table(MMM)

dir.create('classifier', showWarnings = F)
data_unknown = MMM[MMM$Lifestyle %in% c('Unknown')]
data_unknown$Lifestyle = NULL
write.table(data_unknown, 'classifier/data2classify', row.names = T, quote = F, sep = '\t')
data_known = MMM[MMM$Lifestyle %!in% c('Unknown')]




###iterate trough lifestyles
unique_lifestyles = unique(data_known$Lifestyle)
for (i in unique_lifestyles){
  count = length(data_known[data_known$Lifestyle %in% i]$Lifestyle)
  if (count < 9 ){
    next
  }else{
    mtrx = data_known
    mtrx[mtrx$Lifestyle %!in% i]$Lifestyle = '0'
    mtrx[mtrx$Lifestyle %in% i]$Lifestyle = '1'
    filename = paste0('classifier/data_classifier_',i)
    write.table(mtrx, filename, row.names = T, quote = F, sep = '\t')
    cmnd = paste0('python classifier_src/FS_RF2.py -i ', filename, ' -u classifier/data2classify -l ', i)
    system(cmnd)
  }
}


###Create mapping_file augmented
new_mapping_file = data.frame(Sample = character(), Lifestyle = character())
predicted_files = list.files('classifier/', pattern = 'predicted_')

for (i in predicted_files){
  file = data.table(read.csv(paste0('classifier/',i), row.names = 1))
  file = file[file$X1 >0.8]
  lifestyle = str_remove(i, 'predicted_')
  lifestyle = str_remove(lifestyle, '.csv')
  file$Lifestyle = lifestyle
  file = file[,c('Sample', 'Lifestyle')]
  new_mapping_file = rbind(new_mapping_file, file)
}

new_mapping_file = data.table(merge(mapping_file, new_mapping_file, by = 'Sample', all.x = T))
new_mapping_file[new_mapping_file$Lifestyle.x %in% 'Unknown']$Lifestyle.x = new_mapping_file[new_mapping_file$Lifestyle.x %in% 'Unknown']$Lifestyle.y
new_mapping_file = new_mapping_file[,1:2]
colnames(new_mapping_file)[2] = 'Lifestyle'
new_mapping_file[is.na(new_mapping_file$Lifestyle)]$Lifestyle = 'Unknown'
write.table(new_mapping_file, 'mapping_file_augmented.txt', row.names = F)
print('New mapping_file saved as: mapping_file_augmented.txt')

