library(stringr)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
#args[1] is the id2tags input
#args[2] is the prokka_headers.txt input
#args[3] is the output




id2tags<- read.delim2(args[1], header = F)
id2tags<- id2tags[,1:2]
colnames(id2tags)<- c('geneid', 'gene')

prokka_header <- read.delim(args[2], header = F)
prokka_header$prokkadescription <- prokka_header$V1
prokka_header$remove <- sapply(strsplit(prokka_header$V1," "), `[`, 1)

prokka_header$prokkadescription <- str_remove(prokka_header$prokkadescription, prokka_header$remove)
prokka_header$prokkadescription <- sub('.', '', prokka_header$prokkadescription)

prokka_header$remove <- sub('.', '', prokka_header$remove)

prokka_header <- prokka_header[,c(3,2)]
colnames(prokka_header)[1] <- 'gene'

##Merge id2tags with prokka headers
prokka_annotation <- merge(id2tags, prokka_header, by= 'gene', all.x = T)
prokka_annotation <- prokka_annotation[,2:3]
prokka_annotation<- prokka_annotation %>% distinct(geneid, .keep_all = TRUE)

write.table(prokka_annotation, args[3], row.names = F)
