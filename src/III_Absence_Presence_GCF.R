############################
### John Baker Hernandez ###
### johnbakerh@gmail.com ###
###        2020          ###
###  Leiden University   ###
### Bioinformatics Group ###
############################
args = commandArgs(trailingOnly=TRUE)
'R Script for Creation of BGC absence presence per genome from Absence_Presence Table from BiG-SCAPE'
# -----------------------------------------------------------------------------
# Absence Presence Matrix Creation
# Abs_Pres Table Created with Python Script ()
# -----------------------------------------------------------------------------

#Creation of Co_Ocurrence Matrix
'Clustering Matrix With BGC from MIBIG'
abs_pres_matrix <- read.delim(args[1], header= TRUE, sep=",")
abs_pres_matrix <-abs_pres_matrix[, c(2,3)] #removal of first column

#Modify GCF.No and Genome Columns
abs_pres_matrix <- abs_pres_matrix[!grepl("BGC", abs_pres_matrix$Genome),]
abs_pres_matrix$Genome<- gsub("BGC", "_BGC", abs_pres_matrix$Genome)
abs_pres_matrix$GCF.No<- gsub("[a-zA-Z ]", "", abs_pres_matrix$GCF.No) #Removes BiG-SCAPE Class from GCF Column
abs_pres_matrix$GCF.No<- gsub("-_", "", abs_pres_matrix$GCF.No)
abs_pres_matrix$GCF.No <- sub("^", "GCF", abs_pres_matrix$GCF.No) #Appends GCF to the Front of the Number

#Separating Columns for Renaming
abs_pres_matrix$genome2 <- sapply(strsplit(as.character(abs_pres_matrix$Genome),"_"), `[`, 1)
abs_pres_matrix$genome3 <- sapply(strsplit(as.character(abs_pres_matrix$Genome),"_"), `[`, 2)
abs_pres_matrix$genome4 <- sapply(strsplit(as.character(abs_pres_matrix$Genome),"_"), `[`, 3)
abs_pres_matrix$BGC.Region <- sapply(strsplit(as.character(abs_pres_matrix$Genome),"_"), `[`, 4)

#Merging Previous Columns
abs_pres_matrix$Genome <- paste0(abs_pres_matrix$genome2, "_",
                                 abs_pres_matrix$genome3, "_",
                                 abs_pres_matrix$genome4)

#Removal of Final Columns
abs_pres_matrix<- abs_pres_matrix[, -c(3,4,5,6)]

binary_table <- table(abs_pres_matrix)
binary_table[binary_table>0] <- 1
binary_table<- as.data.frame.matrix(binary_table)

columns <- colnames(binary_table)
new_columns <- paste0(sapply(strsplit(columns,"_"), `[`, 1),'_', sapply(strsplit(columns,"_"), `[`, 2), '_', sapply(strsplit(columns,"_"), `[`, 3))
colnames(binary_table) <- new_columns


write.table(binary_table,args[2])
