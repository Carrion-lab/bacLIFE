####################################
###     Guilermo Guerrero        ###
### G.GuerreroEgido@nioo.knaw.nl ###
###           2021               ###
###     Leiden University        ###
###                              ###
####################################


###This script contains several functions used in the shiny app

`%!in%` = Negate(`%in%`)

##Download fasta sequences
extract_sequences_fasta <- function(matrix, all_proteins, db, bacteria, id, names_equivalence){
  if (db == 'COG'){
    extract_genes <- matrix[matrix$cogid %in% id,]
  }
  if (db == 'KEGG'){
    extract_genes <- matrix[matrix$keggid %in% id,]
  }
  if (db == 'PFAM'){
    extract_genes <- matrix[matrix$pfamid %in% id,]
  }
  if (db == 'DBCAN'){
    extract_genes <- matrix[matrix$dbcanid %in% id,]
  }
  if (db == 'MCL'){
    extract_genes <- matrix[matrix$clusters %in% id,]
  }
  
  #names_equivalence = read.table('names_equivalence.txt', header = T)
  names_equivalence$MicroLife_name = str_remove(names_equivalence$MicroLife_name, '_O.fna')
  names_equivalence$Full_name = str_remove(names_equivalence$Full_name, '_O.fna')

  gene_list = character()
  for (i in 1:nrow(extract_genes)){
    gene_list = paste0(gene_list,extract_genes$descriptions[i] )
  }
  
  gene_list = strsplit(gene_list, ',')[[1]]
  
  if (bacteria %in% 'All'){
    to_extract = all_proteins[gene_list]
  }
  if (bacteria %!in% 'All'){
    df = data.frame(fasta_name = gene_list, bacteria_name = sapply(strsplit(gene_list,"[|]"), `[`, 1))
    df = merge(df, names_equivalence, by.x = 'bacteria_name', by.y = 'MicroLife_name')
    df_sub = data.table(df)
    gene_list_one_bacteria = df_sub[df_sub$Full_name %in% bacteria]$fasta_name
    to_extract = all_proteins[gene_list_one_bacteria]
  }
  
  dir.create('fasta_sequences')
  to_extract = to_extract[names(to_extract) %in% NA == FALSE]
  filename = paste0('fasta_sequences/fasta_extract_', id, '_', bacteria, '.fasta')
  ape::write.FASTA(to_extract, filename)
  
  return(filename)
}



#Remove rows with all colums same value (core genes)
remove_constant_columns <- function(matrix){
  n_samples <- ncol(matrix)
  keep <- apply(matrix[2:n_samples], 1, function(x) length(unique(x[!is.na(x)])) != 1)
  new_matrix <-  matrix[keep,]
  return(new_matrix)
}

##MAKE PCOA 
pcoa_function <- function(matrix){
  d <- parallelDist(t(matrix), method = "dice", threads = 30)
  k <- 2
  d[is.na(d)] <- 0
  pcoa <- cmdscale(d, k=k, eig=T)
  eig <- as.numeric(pcoa$eig)
  eig <- 100 * (eig/sum(eig))
  pcoa_data <- data.frame(Sample=rownames(pcoa$points), X= pcoa$points[,1], Y=pcoa$points[,2])
  return_list <- list(pcoa_data, eig)
  
  return(return_list)
}

#Make PCA
pca_function <- function(matrix){
  #matrix <- scale(t(matrix), scale = T, center = T)
  pca<- prcomp(t(matrix), center= T)
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100,1)
  pca.data <- data.frame(Sample=rownames(pca$x), X= pca$x[,1], Y=pca$x[,2])
  return_list <- list(pca.data, pca.var.per)
  #pca_output <- list( pca,pca.var,pca.var.per,pca.data)
  return(return_list)
}

#MAKE BARPLOT
barplot_function <- function(pca.input, db){
  png(paste0('plots/PCs_barplots/', db,   "_barplot_PCs.png")) 
  barplot(pca.input[[3]], main= paste0(db, 'PCs variance'), xlab = 'Principal component', ylab = 'Percent Variation')
  dev.off()
  
  
  return(plot)
}

###Extract core ids by choice (depending on the matrix you give it)
core_groups <- function(matrix, input){
  colnames(matrix)[1] <- 'id'
  vect<- c()
  for (i in names(input)){
    add <- input[[i]]
    vect <- c(vect, add)
    }
  matrix <- matrix[matrix$id %in% vect,]
  return(matrix)
}

###Extract cog processes from hierarchy
cog_options <- function(cog_agg_matrix, cog_extended_annotations){
  mer <- merge(cog_agg_matrix[,1:2], cog_extended_annotations, by = 'id')
  processes <- unique(mer$processes)
  return(processes)
}



###Get COG hierarchy of core genes
all_core_genes_cog <- function(core_genes, cog_extended_annotations){
  
  core_genes_cog <- core_genes[,c(2:n_samples, (n_samples + 6))]
  core_genes_cog <- data.table(core_genes_cog)
  core_genes_cog[is.na(core_genes_cog$cogid)]$cogid <- 'Unknown'
  core_genes_cog <- as.data.frame(core_genes_cog)
  
  core_genes_cog <- data.frame(id = core_genes_cog$cogid, counts = 1) 
  
  core_genes_cog_agg <- aggregate(core_genes_cog[,2], by = list(id = core_genes_cog[,1]), FUN=sum)
  
  core_hierarchy <- merge(core_genes_cog_agg, cog_extended_annotations, by = 'id')
  
  core_hierarchy_agg <- aggregate(core_hierarchy[,2], by = list(processes = core_hierarchy[,3]), FUN = sum)
  colnames(core_hierarchy_agg)[1] <- 'COG_processes'
  
  core_hierarchy_agg$x <- core_hierarchy_agg$x / sum(core_hierarchy_agg$x)
  return(core_hierarchy_agg)
  
}


###Get cog hierarchy of all cd-hit clusters
cog_groups <- function(cog_agg_matrix, cog_extended_annotations, mapping_file, column){
  y_limit <- max(colSums(cog_agg_matrix[,2:ncol(cog_agg_matrix)]))
  M <- merge(cog_agg_matrix, cog_extended_annotations, by = 'id')
  n_samples <- grep("general", colnames(M)) -2
  M_agg <- aggregate(M[,2:(n_samples)], by = list(id = M[,'processes']), FUN = sum)
  mm <- reshape2::melt(M_agg, id='id')
  mmm <- merge(mm, mapping_file, by.x = 'variable', by.y = 'Sample')
  agg <- aggregate(mmm[,3], by = list(variable = mmm$variable), FUN= sum)
  colnames(agg)[2] <- 'n_genes'
  mmm_agg <- merge(mmm, agg, by = 'variable')
  mmm_agg$percent <- (mmm_agg$value / mmm_agg$n_genes) * 100
  cols <- calculate_colors(sort(unique(mmm_agg$id)), df_colors)
  
  plot <- ggplot(mmm_agg[order(mmm_agg[,column],decreasing=TRUE),], aes(fill=id, y=value, x=variable)) + 
    geom_bar( stat="identity")+ theme_bw() + 
    theme(legend.text = element_text(size = 7),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="bottom") + 
    facet_wrap(as.formula(paste("~", column)),  scales = "free") + ylab('Counts') +
    scale_y_continuous(limits=c(0,y_limit)) +  scale_fill_manual(values = cols) + theme(strip.background = element_rect(colour = 'black', fill = 'white'))
  
  return(plot)
}

###Get cog hierarchy of all cd-hit clusters in % values
cog_groups_percent <- function(cog_agg_matrix, cog_extended_annotations, mapping_file, column, process){
  y_limit <- max(colSums(cog_agg_matrix[,2:ncol(cog_agg_matrix)]))
  M <- merge(cog_agg_matrix, cog_extended_annotations, by = 'id')
  n_samples <- grep("general", colnames(M)) -2
  M_agg <- aggregate(M[,2:(n_samples)], by = list(id = M[,'processes']), FUN = sum)
  mm <- reshape2::melt(M_agg, id='id')
  mmm <- merge(mm, mapping_file, by.x = 'variable', by.y = 'Sample')
  agg <- aggregate(mmm[,3], by = list(variable = mmm$variable), FUN= sum)
  colnames(agg)[2] <- 'n_genes'
  mmm_agg <- merge(mmm, agg, by = 'variable')
  mmm_agg$percent <- (mmm_agg$value / mmm_agg$n_genes) * 100
  cols <- calculate_colors(sort(unique(mmm_agg$id)), df_colors)
  
  plot_1 <- ggplot(mmm_agg[order(mmm_agg[,column],decreasing=TRUE),], aes(fill=id, y=value, x=variable)) + 
    geom_bar( stat="identity") + theme_bw() + 
    theme(legend.text = element_text(size = 7),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="bottom") + 
    facet_wrap(as.formula(paste("~", column)),  scales = "free") + ylab('Counts') +
    scale_y_continuous(limits=c(0,y_limit)) +  scale_fill_manual(values = cols) + theme(strip.background = element_rect(colour = 'black', fill = 'white'))
  
  if (process == 'All'){
    to_plot <- mmm_agg
    cols <- calculate_colors(sort(unique(to_plot$id)), df_colors)
    plot_2 <- ggplot(to_plot[order(to_plot[,column],decreasing=TRUE),], aes(fill=id,  x=variable)) + 
      geom_bar(aes(y=percent), stat="identity") + geom_line(aes(y=n_genes/100,group=1), size = 1., color= 'black')+ theme_bw()+ 
      theme(legend.text = element_text(size = 7),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="bottom") + 
      facet_wrap(as.formula(paste("~", column)),  scales = "free") + 
      scale_y_continuous(limits=c(0,100), sec.axis = sec_axis(~.*100, name="number of genes")) +  scale_fill_manual(values = cols) + theme(strip.background = element_rect(colour = 'black', fill = 'white'))
  }else{
    to_plot <- mmm_agg[mmm_agg$id == process,]
    y_max <- max(to_plot$percent)
    plot_2 <- ggplot(to_plot[order(to_plot[,column],decreasing=TRUE),], aes(fill=id, y=percent, x=variable)) + 
      geom_bar( stat="identity")+ theme_bw() + 
      theme(legend.text = element_text(size = 7),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="bottom") + 
      facet_wrap(as.formula(paste("~", column)),  scales = "free") + 
      scale_y_continuous(limits=c(0,y_max)) +  scale_fill_manual(values = c('black')) + theme(strip.background = element_rect(colour = 'black', fill = 'white'))
    
  }
  return_list <- list(plot_1, plot_2)
  return(return_list)
}


cog_core_groups <- function(matrix, input, cog_extended_annotations, process){
  lista <- list()
  for ( i in names(input)){
    df <- data.frame(id= input[[i]], count = 1 )
    cog <- matrix[,c(1, (n_samples + 6))]
    colnames(cog)[1]<- 'id'
    merg <- merge(df, cog, by = 'id')
    merg <- data.table(merg)
    merg[is.na(merg$cogid)]$cogid <- 'Unknown'
    merg <- as.data.frame(merg)
    
    
    merg <- aggregate(merg[,2], by= list(merg$cogid), FUN= sum)
    colnames(merg) <- c('id', 'counts')
    merg <- merge(merg, cog_extended_annotations, by.x = 'id',)
    merg <- aggregate(merg[,2], by= list(processes = merg$processes), FUN = sum)
    merg$x <- (merg$x/sum(merg$x))*100
    colnames(merg)[2] <- i
    lista[[i]] <- merg
  }
  
  df_general <- lista[[1]]
  for(i in 2:length(lista)){
    df_general <- merge(df_general, lista[[i]], all=T)
  }
  
  df_melted <- reshape2::melt(df_general, id = 'processes' )
  
  if (process== 'All'){
    cols <- calculate_colors(sort(unique(df_melted$processes)), df_colors)
    plot <- ggplot(df_melted, aes(x=variable, y = value, fill = processes)) + geom_bar( stat = 'identity') + theme_bw()  + 
      theme(axis.text.x = element_text( angle = 45, hjust =1) ,legend.position="bottom", legend.text = element_text(size = 7)) + ylab('percent %') +
      scale_fill_manual(values = cols)
  }else{
    to_plot <- df_melted[df_melted$processes == process,]
    plot <- ggplot(to_plot, aes(x=variable, y = value, fill = processes)) + geom_bar( stat = 'identity') + theme_bw() + 
      theme(axis.text.x = element_text( angle = 45, hjust =1) ,legend.position="bottom", legend.text = element_text(size = 7))+ ylab('percent %')+
      scale_fill_manual(values = c('black'))
  }
  
  return(plot)
}


#MAKE PCA PLOT
pca_plot_function <- function(pca.input, annotation, groups){
  to_plot <- merge(pca.input[[4]], groups, by = 'Sample')
  plot <- ggplot(to_plot, aes(x=X,y=Y,color=Lifestyle)) + 
    geom_point() + xlab(paste0('PC1 - ', pca.input[[3]][1], '%', sep='')) + 
    ylab(paste0('PC2 - ', pca.input[[3]][2], '%', sep=''))+
    ggtitle(paste0('PCA ', annotation)) + theme_bw()
  return(plot)
}

#MAKE LOADINGS PCA PLOT (not used in app)
pca_loadings_plot_function <- function(pca.input, annotation){
  to_plot <- pca.input[[1]]
  to_plot <- as.data.frame(to_plot$rotation)
  to_plot$Name <- rownames(to_plot)
  plot <- ggplot(to_plot, aes(x= `PC1` ,y=`PC2`, label=Name)) + 
    geom_point() + xlab(paste0('PC1 - ', pca.input[[3]][1], '%', sep='')) + 
    ylab(paste0('PC2 - ', pca.input[[3]][2], '%', sep=''))+ theme_bw() +
    ggtitle(paste0('PCA variable importance ', annotation)) +
    geom_text(aes(label=ifelse(PC1> tail(head(sort(to_plot$PC1,decreasing=TRUE), n = 10), n=1) ,as.character(Name),'')),hjust=0,vjust=0) +
    geom_text(aes(label=ifelse(PC2 > tail(head(sort(to_plot$PC2,decreasing=TRUE), n = 10), n=1) ,as.character(Name),'')),hjust=0,vjust=0)
  return(plot)
}



aggregate_matrix <- function(matrix, db, groups, n_samples){
  
  matrix = aggregate(matrix[,1:(n_samples - 1)], by = list(id = matrix[,n_samples]), FUN=sum)
  matrix<- remove_constant_columns(matrix)
  
  rownames(matrix) <- matrix[,1] 
  matrix[,1] <- NULL
  matrix <- as.matrix(matrix) 
  return(matrix)
}

##Call functions for plots
create_pca <- function(matrix){
  
  pca <- pca_function(matrix)
  
  return(pca)
  
}

###Create upset plot (not used in app)
upset_plot <- function(matrix, db){
  rownames(matrix) <- matrix[,1]
  matrix[,1]<- NULL
  samples <- ncol(matrix)
  matrix[2:samples] <- apply(matrix[2:samples], 2, function(x) ifelse(x > 1, 1, x))
  my_list<- list()
  for (i in 2:samples){
    sample<- colnames(matrix)[i]
    ids_present <- rownames(matrix[matrix[,sample] == 1,])
    my_list[[sample]] <- ids_present
  }
  pdf( file = paste0('corepan_analysis/upset_plots/' ,db,'_upset_plot.pdf'), height = 50 , width = 35)
  plot<- upset(fromList(my_list), order.by = "freq", nsets = (samples))
  print(plot)
  dev.off()
  return(plot)
}

###create upset data object
general_upset_plot <- function(matrix, groups,db, column){
  colnames(matrix)[1] <- 'id'
  lifestyles <- unique(groups[, column])
  my_list <- list()
  for (i in lifestyles){
    samples <- groups[groups[,column] == i,]  
    samples_names <- samples[,1]
    match <-samples_names[samples_names %in% colnames(matrix)]
    subset<- matrix[,c('id', match)]
    n_samples <- ncol(subset)
    subset_core_genes <- core_genes_extraction(subset, n_samples)
    subset_core_genes <- subset_core_genes[[1]]
    
    ids_present <- subset_core_genes[,'id']
    my_list[[i]] <- ids_present
    
  }
  #pdf(file = paste0('corepan_analysis/upset_plots/' ,db,'_general_upset_plot.pdf'), paper = 'USr')
  #plot <- upset(fromList(my_list), order.by = "freq")
  #print(plot)
  #dev.off()
  return(my_list)
}

####Function extract core genes
core_genes_extraction <- function(matrix, n_samples){
  zeros <- rowSums(matrix[2:n_samples] == 0)
  matrix$n_zeros <- zeros
  core_genes <- matrix[matrix$n_zeros < (n_samples - (n_samples*0.9)),]
  non_core_genes <- matrix[matrix$n_zeros > (n_samples - (n_samples*0.9)),]
  core_genes <- core_genes[1:(length(core_genes)-1)]
  return_list <- list(core_genes, non_core_genes)
  return(return_list)
}


####Function extract singletons genes
extract_singletons <- function(matrix){
  n_samples <- grep("completeness", colnames(matrix)) - 1
  zeros <- rowSums(matrix[2:n_samples] == 0)
  matrix$n_zeros <- zeros
  singletons <- matrix[matrix$n_zeros == (n_samples - 2),]
  singletons <- singletons[1:(length(singletons)-1)]
  no_singletons <- matrix[!(matrix$n_zeros == (n_samples - 2)),]
  no_singletons <- no_singletons[1:(length(no_singletons)-1)]
  
  return_list <- list(singletons, no_singletons)
  return(return_list)
  
}


####Function extract genes with more than 90 % of 0s
remove_genes_with_0s <- function(matrix){
  n_samples <- grep("completeness", colnames(matrix)) - 1
  zeros <- rowSums(matrix[2:n_samples] == 0)
  matrix$n_zeros <- zeros
  no_singletons <- matrix[matrix$n_zeros < (n_samples - (n_samples*0.01)),]
  no_singletons <- no_singletons[1:(length(no_singletons)-1)]
  
  return(no_singletons)
  
}


####Function extract singletons per sample and writet them down in corepan_anaysis folder

singleton_per_sample <- function(matrix){
  n_samples <- grep("completeness", colnames(matrix)) - 1
  samples <- colnames(matrix[2:n_samples])
  
  dir.create('corepan_analysis/singletons_per_sample')
  names_list<- c()
  n_singletons_list <- c()
  for (i in 1:length(samples)){
    sample_column <- matrix[,samples[i]]
    keep <- sample_column != 0
    individual_singletons <- matrix[keep, c(1, grep(samples[i], colnames(matrix)), (n_samples + 1):(n_samples + 12))]
    to_print <- individual_singletons[,c(1,5,2,6:14)]
    colnames(to_print)[3] <- 'gene_counts'
    names_list<- c(names_list, samples[i])
    n_singletons_list<- c( n_singletons_list, nrow(to_print))
    write.table(to_print, paste0('corepan_analysis/singletons_per_sample/',samples[i], '_singletons.txt'), row.names = F)
  }
  df <- data.frame(Sample = names_list, singletons = n_singletons_list)
  df <- df[order(-df$singletons),]
  return(df)
  
}



###CLUSTERING###
###          ###

#HDBSCAN
clustering_hdbscan <- function(pca_data, min){
  cl <- hdbscan(pca_data[,2:3],min)
  pca_data$cluster <- factor(cl$cluster)
  pca_data <- data.table(pca_data)
  pca_data[pca_data$cluster == 0]$cluster <- 'Outliers'
  pca_data <- as.data.frame(pca_data)
  #pca_data <- merge(pca_data, mapping, by.x = 'Sample', by.y = 'Sample')
  return(pca_data)
  
}

#KNN
clustering_knn<- function(pca_data,  k){
  res.km <- kmeans(pca_data[,2:3], k, nstart = 25)
  pca_data$cluster <- factor(res.km$cluster)
  #pca_data <- merge(pca_data, mapping, by.x = 'Sample', by.y = 'Sample')
  return(pca_data)
}

##Extract options from mapping file
extract_options <- function(groups, column){
  options <- unique(groups[,column])
  return(options)
}
###Assign colors to a vector (really useful)
assign_colors <- function(vector){
  unique_values <- unique(vector)
  unique_colors <- glasbey(n= length(unique_values))
  df <- data.table(data.frame(vector = vector))
  df$colors <- 'empty'
  
  
  for (i in 1:length(unique_values)){
    df[df$vector %in% unique_values[i]]$colors <- unique_colors[i]
  }
  return(df$colors)
}

##create hierarchical clustering
dendogram <- function(mtrx){
  hc<- hclust(parallelDist(t(mtrx), method = "dice", threads = 30))
  return(hc)
}

### match names of matrix with mapping file
match_sample_names <- function(matrix,mapping_file, n_samples){
  colnames_mapping <- mapping_file$Sample
  colnames_matrix <- colnames(matrix)[2:n_samples]
  
  for(i in colnames_mapping){
    x <- stringdist(i, colnames_matrix)
    position <- x %>% match(x = min(x))
    colnames_matrix[position] <- i
  }
  
  
  return(colnames_matrix)
  
}


calculate_colors <- function(vector, df_colors){
  df <- data.frame(processes = vector)
  M <- merge(df,df_colors, by = 'processes')
  return(M$color)
}

####Extract sequences of the id chosen, i.e id = 'COG4669'
extract_sequences <- function(matrix, all_proteins, db, id){
  
  if (db == 'COG'){
    extract_genes <- matrix[matrix$cogid %in% id,]
  }
  if (db == 'KEGG'){
    extract_genes <- matrix[matrix$keggid %in% id,]
  }
  if (db == 'PFAM'){
    extract_genes <- matrix[matrix$pfamid %in% id,]
  }
  if (db == 'DBCAN'){
    extract_genes <- matrix[matrix$dbcanid %in% id,]
  }
  
  
  genes_id <- extract_genes[,c(1, grep("gene", colnames(extract_genes)))]
  genes_id <- separate_rows(genes_id, gene, sep = ';')
  
  seqs <- all_proteins[genes_id$gene]
  dir.create('MSA')
  write.FASTA(seqs, 'MSA/selected_sequences.fasta')
  
  
}


###Perform MSA alignment
msa_alignment <- function(){
  mySequences <- readAAStringSet('MSA/selected_sequences.fasta')
  myFirstAlignment <- msa(mySequences)
  msaPrettyPrint(myFirstAlignment, output="asis", alFile = 'MSA/alignment.fasta', showNames="none",
                 showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
  
}

##Create tree from fasta MSA file
tree_from_msa <- function(){
  alignment <- readAAStringSet('MSA/alignment.fasta')
  
  ###Plot tree
  d <- as.dist(stringDist(alignment, method = "hamming")/width(alignment)[1])
  tree <- bionj(d)
  p <- ggtree(tree) + geom_tiplab(size = 1.9)
  return_list <- list(tree, p)
  return(return_list)
}


###Load tree into cytoscape
network_cytoscape <- function(tree, mapping_file, column, id, db){
  ig <- as.igraph.phylo(tree)
  ig <- set_edge_attr(ig,'distance', value=tree$edge.length) # set distances as edge attributes
  
  
  ###Get nodes from igraph
  nodes <- get.data.frame(ig, what= c("vertices" ))
  nodes2<- nodes
  
  
  ###Add node attributes
  
  
  #Extract genes of all gene clusters associated with that COG
  
  if (db == 'COG'){
    extract_genes <- matrix[matrix$cogid %in% id,]
  }
  if (db == 'KEGG'){
    extract_genes <- matrix[matrix$keggid %in% id,]
  }
  if (db == 'PFAM'){
    extract_genes <- matrix[matrix$pfamid %in% id,]
  }
  if (db == 'DBCAN'){
    extract_genes <- matrix[matrix$dbcanid %in% id,]
  }
  
  
  genes_id <- extract_genes[,c(1, grep("descriptions", colnames(extract_genes)))]
  
  #Get counts of each cluster
  store <- c()
  for (i in 1:nrow(genes_id)){
    counts <- length(str_split(genes_id$descriptions, ';')[[i]])
    store <- c(store, counts)
  }
  genes_id$n_genes <- store
  
  #Get cluster proportions
  genes_id_2 <- separate_rows(genes_id, descriptions, sep = ';')
  genes_id_2$n_genes <- NULL
  genes_id_2$descriptions <- sapply(strsplit(genes_id_2$descriptions,"[|]"), `[`, 1)
  genes_id_2$descriptions <- str_replace_all(genes_id_2$descriptions, '-', '.')
  genes_id_2_attributes <- merge(genes_id_2, mapping_file, by.x = 'descriptions', by.y = 'Sample', all.x = T)
  
  #column <- 'clusters_knn'
  cluster_rates <- reshape2::dcast(genes_id_2_attributes, as.formula(paste('clusters ~ ', column)) )
  
  
  ####Merge cluster rates with cluster counts and get them ready to load in cytoscape
  genes_id$descriptions <- NULL
  M <- merge(genes_id, cluster_rates, by = 'clusters')
  
  genes_id <- extract_genes[,c(1, grep("gene", colnames(extract_genes)))]
  M2 <- merge(M, genes_id, by = 'clusters')
  
  M3 <- merge( nodes2,M2, by.x = 'name', by.y = 'gene', all.x = T)
  
  
  row.names(M3) <- M3$name
  M3$name <- NULL
  M3$n_genes[is.na(M3$n_genes)] <- 0
  M3$n_genes <- M3$n_genes/sum(M3$n_genes)
  ###Load network in Cytoscape
  
  createNetworkFromIgraph(ig, title="phylotree", collection = "phylotree")
  #setEdgeLabelMapping('distance')
  
  layoutNetwork(paste('force-directed',
                      'edgeAttribute="distance"',
                      'type="1 - normalized value"',
                      'defaultSpringCoefficient=5E-4',
                      'defaultSpringLength=50',
                      sep = ' '))
  
  loadTableData(M3)
  
  setNodeCustomPieChart(colnames(M3)[3:ncol(M3)])
  
  lockNodeDimensions(TRUE)
  setNodeShapeDefault('ELLIPSE')
  setNodeSizeMapping('n_genes', c(0.01,1), c(20,200))
  
  
  
}


