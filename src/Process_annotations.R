####################################
###     Guilermo Guerrero        ###
### G.GuerreroEgido@nioo.knaw.nl ###
###           2021               ###
###     Leiden University        ###
###                              ###
####################################

###This script merge annotations with the cd-hit abs/pres matrix and produces the MEGAMATRIX


args = commandArgs(trailingOnly=TRUE)
##arg[1] is the matrix from cdhit
##arg[2] is the cog annotatuion
##arg[3] is the kegg annotatuion
##arg[4] is the hmm annotation
##arg[5] is the clusterVSid file
##arg[6] is the prokka annotations
##arg[7] is the output name of matrix with all annotations





'Reading input files'
matrix <- read.table(args[1], header = T)
clustervsgene <- matrix[,c('clusters', 'gene')]

###COG annotations-
cog <- read.table(args[2], header = T)

####KEGG annotations
kegg <- read.table(args[3], header = T)

###Prokka annotations
prokka <- read.table(args[5], header = T)

'Merge annotations COG/KEGG'
###Merge annotations
all_annotations <- merge(kegg, cog, by= 'gene', all = T)

'Merge with HMM annotations'
##Merge with hmm annotations
hmm_annotations <- read.table(args[4], header = T)
all_annotations <- merge(all_annotations, hmm_annotations, by= 'gene', all = T)

'Merge with prokka annotation'
##Merge with prokka
all_annotations <- merge(all_annotations, prokka, by.x= 'gene', by.y = 'geneid', all = T)



###Merge with absence_presence_matrix

'Merge with absence_presence_matrix'
all_annotations <- merge(clustervsgene, all_annotations, by = 'gene', all=T)
all_annotations$gene<- NULL
column_number= ncol(matrix)
matrix_info = matrix[,c(column_number -1 , column_number )]
matrix <- matrix[1:(column_number-2)]
matrix$completeness <- 'True'
matrix <- cbind(matrix, matrix_info)


MEGAMATRIX <- merge(matrix, all_annotations, by= 'clusters', all.x=T)
'Write MEGAMATRIX'
write.table(MEGAMATRIX, args[6], row.names = F)

