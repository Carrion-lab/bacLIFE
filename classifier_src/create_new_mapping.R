###Join mapping file with the classification of the RF

library(stringr)
library(dplyr)
library(readr)
library(data.table)
source('src/corepan_functions.R')
source('src/stats_functions.R')
source('src/text_shiny.R')


`%!in%` = Negate(`%in%`)




###Load metadata
mapping_file <- read.table('mapping_file_original.txt', header = T)

pp_class = read.table('classifier/predicted_plant_pathogen.csv', header = T, row.names = NULL, sep = ',')
pp_class$X = NULL
pp_class$X0 = NULL
ap_class = read.table('classifier/predicted_opportunistic_pathogen.csv', header = T, row.names = NULL, sep = ',')
ap_class$X = NULL
ap_class$X0 = NULL
env_class = read.table('classifier/predicted_Environmental_Plant_beneficial.csv', header = T, row.names = NULL, sep = ',')
env_class$X = NULL
env_class$X0 = NULL

M = merge(mapping_file, pp_class, by = 'Sample', all.x = T)
colnames(M)[3] = 'plant_pathogen'
M = merge(M, ap_class, by = 'Sample', all.x = T)
colnames(M)[4] = 'opportunistic_pathogen'
M = merge(M, env_class, by = 'Sample', all.x = T)
colnames(M)[5] = 'environmental'

M = data.table(M)
M[M$plant_pathogen > 0.7]$Lifestyle = 'plant_pathogen'
M[M$opportunistic_pathogen > 0.7]$Lifestyle = 'opportunistic_pathogen'
M[M$environmental > 0.7]$Lifestyle = 'Environmental_Plant_beneficial'

write.table(M[,c(1,2)], 'mapping_file_classifier.txt', row.names = F)
