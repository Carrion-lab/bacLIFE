############################
### John Baker Hernandez ###
### johnbakerh@gmail.com ###
###        2020          ###
###  Leiden University   ###
### Bioinformatics Group ###
############################

'R Script for Data Mining and Manipulation of Output Tables from BiG-SCAPE'

'Libraries to be Imported'
library(plyr)
library(dplyr)
library(stats)
library(tidyverse)
library(ggplot2)

#Set Working Directory
#setwd('/Results_I')
#ex:setwd("E:/NACTAR")

args = commandArgs(trailingOnly=TRUE)
"Import Data Tables from BiG-SCAPE"
#Clustering Tables
clustering_df_mix <- read.delim(args[1], header = TRUE) #clustering dataframe

#Network Tables
network_df <- read.delim(args[2], header = TRUE) #general network dataframe
network_df <- network_df [-c(4:11)] #removal of unwanted columns

#Annotation Tables
annotations_df_full<- read.delim(args[3], header = TRUE) #annotations dataframe

'Merge Files for Absence Presence Table'
annotations_merged_df_mix <- merge(annotations_df_full , clustering_df_mix , by.x= "BGC", by.y = "X.BGC.Name") #merge clustering and annotations
annotations_merged_df_mix  <- annotations_merged_df_mix [, c(1,8,6,5)]
#annotations_merged_df_mix  <- annotations_merged_df_mix [, c(1,8,4,5,6)] #removal of unwanted columns

colnames(annotations_merged_df_mix) <- c('BGC Name', 'GCF No', 'Organism', 'BGC Class')

'Extract Rows based on Logical Criteria (User Based Criteria) (Optional by User)'
#sum(network_df_mix$Raw.distance > 0.1)
network_df__distance_filtered <- network_df %>% filter(Raw.distance > 0.10)
network_df_bgc_filtered <- network_df[!grepl("BGC", network_df$Clustername.1),]



#Export Tables to CSV Format
# Write data to txt file: tab separated values
# sep = "\t"
write.table(network_df__distance_filtered, file = args[4], sep = "\t",
            row.names = FALSE, quote = FALSE)
# Write data to csv files:
# decimal point = "." and value separators = comma (",")

write.table(annotations_merged_df_mix, file = args[5], sep = '\t',
            row.names = FALSE, quote = FALSE)




###Obtain GCF annotations list
merged_annotations <- annotations_merged_df_mix
merged_annotations$GCF.No <- paste0( 'GCF',merged_annotations$GCF.No)
merged_annotations<- merged_annotations[,c(2,4)]
merged_annotations <- distinct(merged_annotations)

write.csv(merged_annotations, args[6], row.names = F)