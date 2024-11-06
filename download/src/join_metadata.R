##Script to join metadata

NCBIid = read.table('names_NCBIid.txt', header = F)
NCBIid$V1 <- paste0(NCBIid$V1, '.fa')

clusters = read.table('cluster_membership.txt', header = T)



M <- merge(clusters, NCBIid, by.x = 'Sample', by.y = 'V1')
colnames(M) <- c('NCBIid', 'cluster_membership', 'representative', 'scientific_name')

write.table(M, 'METADATA_MERGED.txt', row.names = F, quote = F)

