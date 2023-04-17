library(stringr)
library(readr)
args = commandArgs(trailingOnly=TRUE)


matrix <- read_delim(args[1], col_names = T, quote = "\"")

matrix <- as.data.frame(matrix)

names_equivalence <- read.table(args[3], header = T)
names_equivalence$Full_name <- str_remove(names_equivalence$Full_name, '_O.fna')
names_equivalence$MicroLife_name <- str_remove(names_equivalence$MicroLife_name, '_O.fna')

n_samples <- grep("completeness", colnames(matrix)) - 1

old_colnames <- data.frame(MicroLife_name = colnames(matrix)[2:n_samples] ) 

M <- merge(old_colnames, names_equivalence, by = 'MicroLife_name')

M <- M[match(M$MicroLife_name,old_colnames$MicroLife_name),]
colnames(matrix)[2:n_samples] <- M$Full_name

write.table(matrix, args[4], row.names = F)


###Rename bigscape table
big_scape_matrix <- read.table(args[2], header = T)


old_colnames <- data.frame(MicroLife_name = colnames(big_scape_matrix) ) 

M <- merge(old_colnames, names_equivalence, by = 'MicroLife_name')

M <- M[match(M$MicroLife_name,old_colnames$MicroLife_name),]
colnames(big_scape_matrix)[2:(ncol(big_scape_matrix) - 1)]  <- M$Full_name

write.table(big_scape_matrix, args[5], row.names = T)

###Create mapping_file
mapping_file = data.frame(Sample = names_equivalence$Full_name, Lifestyle = 'Unknown')
write.table(mapping_file, args[6], row.names = F)
