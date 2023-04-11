library(stringr)
library(dplyr)
library("tibble")

args = commandArgs(trailingOnly=TRUE)
#args[1] is the MEGAMATRIX original input
#args[2] is the new abs/pres matrix
#args[3] is the old mapping file
#args[4] is the new sample lifestyle, default NEW


##LOAD MEGAMATRIX and new abs/pres matrix 
MEGAMATRIX <- read.table(args[1], header = T)
new_matrix <- read.table(args[2], header = T)

###Get colnames
cols_MEGAMATRIX <- colnames(MEGAMATRIX)
cols_new_matrix <- colnames(new_matrix)

##Look which colnames are new
new_col_index <- which(!cols_new_matrix %in% cols_MEGAMATRIX)
new_col <- new_matrix[,c(1,new_col_index)]

##Merge matrices and put new column in specific position
new_MEGAMATRIX <- merge(MEGAMATRIX, new_col, by= 'clusters', all.y = T)

last_column_MEGAMATRIX <- grep("completeness", colnames(MEGAMATRIX)) - 1
new_MEGAMATRIX <- new_MEGAMATRIX %>% relocate(colnames(new_matrix)[new_col_index], .before = colnames(MEGAMATRIX)[last_column_MEGAMATRIX + 1])
last_column_new_MEGAMATRIX <- grep("completeness", colnames(new_MEGAMATRIX)) - 1
new_MEGAMATRIX[, 1:last_column_new_MEGAMATRIX][is.na(new_MEGAMATRIX[, 1:last_column_new_MEGAMATRIX])] <- 0

write.table(new_MEGAMATRIX, 'MEGAMATRIX_noupdated.txt', row.names = F)


###UPDATE MAPPING_FILE

mapping_file <- read.table(args[3], header = T)
new_sample <- colnames(new_matrix)[new_col_index]
new_row <- data.frame(Sample= new_sample, Lifestyle= args[4])

new_mapping_file <- rbind(mapping_file, new_row)

write.table(new_mapping_file, 'corrected_mapping_file.txt', row.names = F)


