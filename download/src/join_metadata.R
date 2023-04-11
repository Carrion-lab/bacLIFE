##Script to join metadata
names_equivalence = read.table('names_equivalence.txt', header = T)

NCBIid = read.table('names_NCBIid.txt', header = F)
NCBIid$V1 <- paste0(NCBIid$V1, '.fa')

clusters = read.table('cluster_membership.txt', header = T)



M <- merge(clusters, NCBIid, by.x = 'Sample', by.y = 'V1')
M <- merge(M, names_equivalence, by.x = 'V2', by.y = 'Full_name', all.x = T)
colnames(M) <- c('scientific_name', 'NCBIid', 'cluster_membership', 'representative', 'MicroLife_name')

write.table(M, 'METADATA_MERGED.txt', row.names = F, quote = F)

##Create mapping_file
mapping_file = data.frame(Sample = M$scientific_name, Lifestyle = 'Unknown')
write.table(mapping_file, 'mapping_file.txt', row.names = F, quote = F)
