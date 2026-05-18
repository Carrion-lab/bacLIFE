library(stringr)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
#args[1] is the id2tags input
#args[2] is the bakta_headers.txt input
#args[3] is the output




id2tags<- read.delim2(args[1], header = F)
id2tags<- id2tags[,1:2]
colnames(id2tags)<- c('geneid', 'gene')

bakta_header <- read.delim(args[2], header = F)
bakta_header$baktadescription <- bakta_header$V1
bakta_header$remove <- sapply(strsplit(bakta_header$V1," "), `[`, 1)

bakta_header$baktadescription <- str_remove(bakta_header$baktadescription, bakta_header$remove)
bakta_header$baktadescription <- sub('.', '', bakta_header$baktadescription)

bakta_header$remove <- sub('.', '', bakta_header$remove)

bakta_header <- bakta_header[,c(3,2)]
colnames(bakta_header)[1] <- 'gene'

##Merge id2tags with bakta headers
bakta_annotation <- merge(id2tags, bakta_header, by= 'gene', all.x = T)
bakta_annotation <- bakta_annotation[,2:3]
bakta_annotation<- bakta_annotation %>% distinct(geneid, .keep_all = TRUE)

write.table(bakta_annotation, args[3], row.names = F)
