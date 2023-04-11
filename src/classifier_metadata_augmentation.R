##############################################
###         Guilermo Guerrero              ###
### g.guerrero.egido@biology.leidenuniv.nl ###
###               2021                     ###
###          Leiden University             ###
###                                        ###
##############################################



##This script takes as input the MEGAMATRIX and the mapping file to augment metadata information with a 80% model confidence
##The output consists in the mapping_file.txt augmented with lifestyle information in genomes labeled as 'Unknown' when confidence is high

library(randomForest)
library(pROC)
library(data.table)
library(readr)
library(parallelDist)
source('src/corepan_functions.R')
source('src/stats_functions.R')
args = commandArgs(trailingOnly=TRUE)
`%!in%` = Negate(`%in%`)




###Load metadata
mapping_file <- read.table(args[1], header = T)

###load and filter matrix
matrix <- read_delim(args[2], col_names = T, quote = "\"")
matrix = as.data.frame(matrix)
n_samples <- grep("completeness", colnames(matrix)) - 1
matrix[,2:n_samples] <- sapply(matrix[,2:n_samples],as.numeric)

#Extract core_genes
dir.create('corepan_analysis')
core_genes_out <- core_genes_extraction(matrix, n_samples)
core_genes <- core_genes_out[[1]]
write.table(core_genes,'corepan_analysis/core_genes.txt', row.names = F)
non_core_genes <- core_genes_out[[2]]

##Extract singletons genes
extract_singleton_genes <- extract_singletons(non_core_genes)
singletons <- extract_singleton_genes[[1]]
no_singletons <- extract_singleton_genes[[2]]

##Remove genes with more than 95% of 0s
clean_matrix <- remove_genes_with_0s(no_singletons)

##Extract singletons for each sample(output automatic in directory singletons_per_sample/ )
singletons_table <- singleton_per_sample(singletons)


##Prepare filtered matrix
n_samples <- grep("completeness", colnames(clean_matrix)) - 1
clean_matrix[2:n_samples] <- apply(clean_matrix[2:n_samples], 2, function(x) ifelse(x > 1, 1, x))
full_matrix <- clean_matrix[,1:n_samples]
full_matrix<- remove_constant_columns(full_matrix)

rownames(full_matrix) <- full_matrix$clusters 
full_matrix$clusters <- NULL
full_matrix <- as.matrix(full_matrix)



#####Apply classifier for metadata augmentation

mapping_file_clean = mapping_file[mapping_file$Lifestyle != 'Unknown',]
mapping_file_unknown = mapping_file[mapping_file$Lifestyle %in% 'Unknown',]

mapping_file_clean <- data.table(mapping_file_clean)
#mapping_file_clean <- mapping_file_clean[-sample(which(Lifestyle=="cepacia_complex"), 286)]

samples_labelled <- mapping_file_clean$Sample


##make Pcoa
d <- parallelDist(t(full_matrix), method = "dice", threads = 30)
k <- 15
d[is.na(d)] <- 0
pcoa <- cmdscale(d, k=k, eig=T)
eig <- as.numeric(pcoa$eig)
eig <- 100 * (eig/sum(eig))

#ncomponents
n_components = sum(eig >= 1)

pcoa_data <- data.frame(Sample=rownames(pcoa$points))
for (i in 1:n_components){
  pcoa_data[,paste0('PC', i)] = pcoa$points[,i]
}

full_pca <- merge(pcoa_data, mapping_file, by = 'Sample')
rownames(full_pca) = full_pca$Sample
full_pca$Sample = NULL
full_pca$Species = NULL

full_pca_to_model = full_pca[samples_labelled,]
full_pca_to_model$Lifestyle <- as.factor(full_pca_to_model$Lifestyle)

full_pca_to_predict = full_pca[mapping_file_unknown$Sample,]

full_pca_to_predict$Lifestyle <- as.factor(full_pca_to_predict$Lifestyle)

###train validation split
set.seed(100)
train <- sample(nrow(full_pca_to_model), 0.7*nrow(full_pca_to_model), replace = FALSE)
TrainSet <- full_pca_to_model[train,]
ValidSet <- full_pca_to_model[-train,]


# Create a Random Forest model with default parameters
model1 <- randomForest(Lifestyle ~ ., data = full_pca_to_model, importance = TRUE, ntree = 500) #sampsize is size of min group *0.8
model1

predictions <- as.data.frame(predict(model1, full_pca_to_predict, type = "prob"))



##Unique lifestyles

unique_lifestyles = unique(mapping_file$Lifestyle)
unique_lifestyles <- unique_lifestyles[!unique_lifestyles %in% c('unknown', 'Unknown')]



###Extract predictions
new_mapping_file = data.table(mapping_file)

for (i in 1:length(unique_lifestyles)){
  idx = grep(unique_lifestyles[i], colnames(predictions))
  pred = rownames(predictions[predictions[,idx] > 0.8,])
  new_mapping_file[new_mapping_file$Sample %in% pred]$Lifestyle = unique_lifestyles[i]
  
}


write.table(new_mapping_file, 'mapping_file_classifier.txt', row.names = F)





