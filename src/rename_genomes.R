library(stringr)
args = commandArgs(trailingOnly=TRUE)


files = list.files(args[1], pattern= '.fna')



df <- data.frame(Full_name = files, Genus = sapply(strsplit(files,"_"), `[`, 1), Species = sapply(strsplit(files,"_"), `[`, 2), Strain = sapply(strsplit(files,"_"), `[`, 3) )

df$Strain <- paste0( 'X',sprintf('%0.5d', 1:nrow(df)))
df$Species <- str_remove(df$Species, '[.]')
df$MicroLife_name <- paste0(df$Genus, '_', str_remove_all(df$Species, '[.]') , '_', df$Strain, '_O.fna')

file.rename(paste0('data/',df$Full_name),paste0( 'data/', df$MicroLife_name))


write.table(df[,c('Full_name', 'MicroLife_name')], 'names_equivalence.txt', row.names = F)
