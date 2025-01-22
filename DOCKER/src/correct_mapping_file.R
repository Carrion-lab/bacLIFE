####################################
###     Guilermo Guerrero        ###
### G.GuerreroEgido@nioo.knaw.nl ###
###           2021               ###
###     Leiden University        ###
###                              ###
####################################



###This script changes the name of the samples to match with the MEGAMATRIX

args = commandArgs(trailingOnly=TRUE)
#args[1] is the mapping_file.txt input
#args[2] is the output (corrected mapping file)

library(stringr)
df <- read.table(args[1], header = T)



df$Sample <- str_remove(df$Sample, '.fna')
df$Sample <- paste0(sapply(strsplit(df$Sample,"_"), `[`, 1), '_', sapply(strsplit(df$Sample,"_"), `[`, 2), '_', sapply(strsplit(df$Sample,"_"), `[`, 3))


write.table(df, args[2], sep = '\t', row.names = F)
