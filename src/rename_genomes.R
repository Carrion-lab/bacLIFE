library(stringr)

message("Renaming genome files and creating names_equivalence.txt")
message("Please note that Capital letters and special characters will be
     removed from the genus and species names.")
message("For example, if the genus name is Bacillus-B,
     it will be renamed to Bacillusb.")

args <- commandArgs(trailingOnly = TRUE)


files <- list.files(args[1], pattern = '.fna')

clean_genus <- function(x) {
  x <- gsub("[^A-Za-z0-9]", "", x)   # keep only letters and numbers
  x <- tolower(x)                    # all lowercase
  x <- paste0(toupper(substr(x, 1, 1)),
              substr(x, 2, nchar(x)))  # capitalize first letter
  return(x)
}
clean_species <- function(x) {
  x <- gsub("[^A-Za-z0-9]", "", x)   # keep only letters and numbers
  x <- tolower(x)                    # all lowercase
  return(x)
}


df <- data.frame(Full_name = files, 
                 Genus = sapply(strsplit(files, "_"), `[`, 1),
                 Species = sapply(strsplit(files, "_"), `[`, 2),
                 Strain = sapply(strsplit(files, "_"), `[`, 3))

df$Genus <- clean_genus(df$Genus)
df$Species <- clean_species(df$Species)

df$baclife_id <- paste0( 'X',sprintf('%0.5d', 1:nrow(df)))

df$bacLIFE_name <- paste0(df$Genus, '_', tolower(str_remove_all(df$Species, '[.]')) , '_', df$baclife_id, '_O.fna')

df$name_reformatted <- paste0(df$Genus,'_',df$Species,'_', df$Strain, '_O.fna')

file.rename(paste0('data/',df$Full_name),paste0( 'data/', df$bacLIFE_name))


write.table(df[,c('Full_name', 'bacLIFE_name')], 'names_equivalence.txt', row.names = F)
