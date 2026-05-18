library(stringr)
library(readr)


args = commandArgs(trailingOnly=TRUE)

message("Renaming MAtrix and bigscape annotation files")

matrix <- read_delim(args[1], col_names = T, quote = "\"", name_repair = "minimal")
matrix <- as.data.frame(matrix)
names_equivalence <- read.table(args[3], header = T)

names_equivalence$Full_name <- str_remove(names_equivalence$Full_name, '_O.fna')
names_equivalence$bacLIFE_name <- str_remove(names_equivalence$bacLIFE_name, '_O.fna')

n_samples <- grep("completeness", colnames(matrix)) - 1

old_colnames <- data.frame(bacLIFE_name = colnames(matrix)[2:n_samples] ) 

M <- merge(old_colnames, names_equivalence, by='bacLIFE_name', all.x=TRUE)
M <- M[match(M$bacLIFE_name,old_colnames$bacLIFE_name),]

# Where Full_name is NA (new genomes not yet in equivalence table), keep bacLIFE_name
M$Full_name <- ifelse(is.na(M$Full_name), M$bacLIFE_name, M$Full_name)

# Report any unmapped genomes
unmapped <- M$bacLIFE_name[M$Full_name == M$bacLIFE_name]
if (length(unmapped) > 0) {
  message("NOTE: ", length(unmapped), " genome(s) not in names_equivalence.txt — kept as-is:")
  message(paste(" ", unmapped, collapse="\n"))
  message("Add them to names_equivalence.txt for full renaming.")
}

colnames(matrix)[2:n_samples] <- M$Full_name

write.table(matrix, args[4], row.names = F)


###Rename bigscape table
big_scape_matrix <- read.table(args[2], header = T, check.names = FALSE)

old_bs_colnames   <- data.frame(bacLIFE_name = colnames(big_scape_matrix))
M2 <- merge(old_bs_colnames, names_equivalence, by='bacLIFE_name', all.x=TRUE)
M2 <- M2[match(old_bs_colnames$bacLIFE_name, M2$bacLIFE_name), ]
M2$Full_name <- ifelse(is.na(M2$Full_name), M2$bacLIFE_name, M2$Full_name)
colnames(big_scape_matrix) <- M2$Full_name
write.table(big_scape_matrix, args[5], row.names=TRUE)

message("Creating mapping file for the Shiny app")

###Create mapping_file
mapping_file = data.frame(Sample = names_equivalence$Full_name, Lifestyle = 'Unknown')
write.table(mapping_file, args[6], row.names = F)
