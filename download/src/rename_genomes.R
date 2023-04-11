library(stringr)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

new_names = read.table('names_NCBIid.txt', header = F)
new_names$V1 = paste0(new_names$V1, '.fa')


files = list.files('genomes_renamed/', pattern= '.fa')



new_names = new_names[new_names$V1 %in% files,]

df <- data.frame(NCBIid =  new_names$V1 ,Full_name = new_names$V2, Genus = sapply(strsplit(new_names$V2,"_"), `[`, 1), Species = sapply(strsplit(new_names$V2,"_"), `[`, 2), Strain = sapply(strsplit(new_names$V2,"_"), `[`, 3) )

df$Strain <- paste0( 'X',sprintf('%0.5d', 1:nrow(df)))

df$MicroLife_name <- paste0(df$Genus, '_', df$Species, '_', df$Strain, '_O.fna')

file.rename(paste0('genomes_renamed/',df$NCBIid),paste0( 'genomes_renamed/', df$MicroLife_name))


write.table(df[,c('Full_name', 'MicroLife_name')], 'names_equivalence.txt', row.names = F)
