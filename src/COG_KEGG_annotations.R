args = commandArgs(trailingOnly=TRUE)
##arg[1] is the annotation file from emapper
##arg[2] is the cog descriptions
##arg[3] is the kegg descriptions
##arg[4] is the cog output
##arg[5] is the kegg output

annotation <- read.delim2(args[1], header = F)
annotation <-annotation[,c(1,5,12)]
colnames(annotation) <- c('id', 'cogid', 'keggid')



#Extract COGs
cog_annotation <- annotation[,c(1,2)]
keep <- grep("COG", cog_annotation$cogid)
cog_annotation <-cog_annotation[keep,]
cog_annotation$cogid <- sapply(strsplit(cog_annotation$cogid,"@"), `[`, 1)
keep <- grep("ar", cog_annotation$cogid, invert = T)
cog_annotation <-cog_annotation[keep,]


##Load COG descriptions
cog_extended_annotations <- read.csv(args[2], header = F)
cog_extended_annotations <- cog_extended_annotations[,c(1,3)]
colnames(cog_extended_annotations)[1] <- 'id'
colnames(cog_extended_annotations)[2] <- 'cog_description'

cog_annotation<- merge(cog_annotation, cog_extended_annotations, by.x = 'cogid', by.y = 'id', all.x = T)
cog_annotation<- cog_annotation[,c(2,1,3)] ##TO WRITE
colnames(cog_annotation)[1]<- 'gene'

write.table(cog_annotation, args[4], row.names = F)


#Extract KEGGs
kegg_annotation <- annotation[,c(1,3)]
keep <- grep("ko", kegg_annotation$keggid)
kegg_annotation <-kegg_annotation[keep,]
kegg_annotation$keggid <- sapply(strsplit(kegg_annotation$keggid,","), `[`, 1)
kegg_annotation$keggid <- sapply(strsplit(kegg_annotation$keggid,":"), `[`, 2)

##KO descriptions
ko_descriptions <- read.delim2(args[3], header = F)
ko_descriptions$V1 <- sapply(strsplit(ko_descriptions$V1,":"), `[`, 2)
colnames(ko_descriptions)<- c('keggid', 'kegg_description')

kegg_annotation <- merge(kegg_annotation, ko_descriptions, by= 'keggid', all.x = T)
kegg_annotation <-kegg_annotation[,c(2,1,3)]
colnames(kegg_annotation)[1]<- 'gene'

write.table(kegg_annotation, args[5], row.names = F)



