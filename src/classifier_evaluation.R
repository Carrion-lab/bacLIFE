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

#Extract singletons genes
extract_singleton_genes <- extract_singletons(non_core_genes)
singletons <- extract_singleton_genes[[1]]
no_singletons <- extract_singleton_genes[[2]]

#Remove genes with more than 95% of 0s
clean_matrix <- remove_genes_with_0s(no_singletons)

#Extract singletons for each sample(output automatic in directory singletons_per_sample/ )
singletons_table <- singleton_per_sample(singletons)


#Prepare filtered matrix
n_samples <- grep("completeness", colnames(clean_matrix)) - 1
clean_matrix[2:n_samples] <- apply(clean_matrix[2:n_samples], 2, function(x) ifelse(x > 1, 1, x))
full_matrix <- clean_matrix[,1:n_samples]
full_matrix<- remove_constant_columns(full_matrix)

rownames(full_matrix) <- full_matrix$clusters 
full_matrix$clusters <- NULL
full_matrix <- as.matrix(full_matrix)


###Generate PCoA
d <- parallelDist(t(full_matrix), method = "dice", threads = 30)
k <- 15
d[is.na(d)] <- 0
pcoa <- cmdscale(d, k=k, eig=T)
eig <- as.numeric(pcoa$eig)
eig <- 100 * (eig/sum(eig))

n_components = sum(eig >= 1)

pcoa_data <- data.frame(Sample=rownames(pcoa$points))
for (j in 1:n_components){
  pcoa_data[,paste0('PC', j)] = pcoa$points[,j]
}


###Iterate the model per each class in the metadata column 'Lifestyle'

unique_lifestyles = unique(mapping_file$Lifestyle)
unique_lifestyles <- unique_lifestyles[!unique_lifestyles %in% c('unknown', 'Unknown')]

prediction_list = list()
AUC_list = list()
ROC_lists = list()


for (i in 1:length(unique_lifestyles)){
  mapping_file_clean = mapping_file[mapping_file$Lifestyle != 'Unknown',]
  mapping_file_unknown = mapping_file[mapping_file$Lifestyle %in% 'Unknown',]
  mapping_file_clean <- data.table(mapping_file_clean)
  mapping_file_clean[mapping_file_clean$Lifestyle %!in% unique_lifestyles[i] ]$Lifestyle = 'Unknown'
  
  

  samples_labelled <- mapping_file_clean$Sample
  full_pca <- merge(pcoa_data, mapping_file_clean, by = 'Sample')
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
  
  # Create a Random Forest model with default parameters and adjust sample size to avoid class imbalance
  if (sum(TrainSet$Lifestyle == unique_lifestyles[i]) < sum(TrainSet$Lifestyle != unique_lifestyles[i])){
    wn = sum(TrainSet$Lifestyle == unique_lifestyles[i]) * 0.8
  }else{
    wn = sum(TrainSet$Lifestyle != unique_lifestyles[i]) * 0.8
  }
  wn = round(wn)
  
  model1 <- randomForest(Lifestyle ~ ., data = TrainSet, importance = TRUE, sampsize = c(wn,  wn))
  model1
  
  predictions <- as.data.frame(predict(model1, ValidSet, type = "prob"))
  
  roc.multi <- multiclass.roc(ValidSet$Lifestyle, as.matrix(predictions), percent= T)
  AUC = auc(roc.multi)
  
  print(paste0('Model evaluation for class: ', unique_lifestyles[i]))
  print(paste0('area under the curve: ',  AUC, '%'))
  
  prediction_list[[i]] = predictions
  AUC_list[[i]] = AUC
  ROC_lists[[i]] = roc.multi
}



###Create and save plots
dir.create('classifier_plots')
for (i in 1:length(unique_lifestyles)){
  AUC_list[[i]]
  ROC = ROC_lists[[i]]
  ROC_animal_pathogen = data.frame(Sensitivity = ROC$rocs[[1]][[1]][["sensitivities"]], Specificity = ROC$rocs[[1]][[1]][["specificities"]])
  png(paste0( 'classifier_plots/' , unique_lifestyles[i], "_classifier_plot.png"), width = 450, height = 400)
  plot(100 - ROC_animal_pathogen$Sensitivity,ROC_animal_pathogen$Specificity,frame = FALSE, type = "l", xlab="1 - Sensitivity", ylab="Specificity", main = paste0(unique_lifestyles[i], ' classifier AUC') , col = '#ff366b' )
  dev.off()
}

'Plots saved in directory classifier_plots/'




