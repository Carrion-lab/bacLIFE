###Script to rename the names
library(stringr)
library(data.table)
sci_names = read.delim2('scientific_names.txt', sep = '\t', header = F)


df = data.frame(NCBI_id = sci_names$V1, sci_names = sci_names$V2, strain_names =  sci_names$V3)
df$Genus = sapply(strsplit(df$sci_names," "), `[`, 1)
df$Species = sapply(strsplit(df$sci_names," "), `[`, 2)
df$strain_names = str_replace_all( df$strain_names, ' ', '.')
df$strain_names = str_replace_all( df$strain_names, '/', '.')
df$strain_names = str_replace_all( df$strain_names, ';', '')
df = data.table(df)
df[df$Species %in% 'sp.']$Species = 'sp'
df$NCBI_id = str_replace_all(df$NCBI_id, '_', '.')

##Write new name
df$new_name = paste0(df$Genus, '_', df$Species,'_', df$strain_names, '|', df$NCBI_id)
df = df[,c(1, 6)]
df$NCBI_id = str_replace(df$NCBI_id, '[.]', '_')

write.table(df, 'names_NCBIid.txt', row.names = F, col.names = F)
