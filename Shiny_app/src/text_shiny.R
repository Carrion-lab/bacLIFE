#####Text for shiny app

introduction <- 'MicroLife Shiny app allows you to interactively explore the gene level and biosynthetic gene level  absence/presence matrices generated in the Snakemake module. An overview of this module is shown below.'

exploration <- 'Principal coordinate analysis based in dice dissimilarity and dendrogram visualization of the gene level and biosynthetic gene level absence/presence matrices annotated with input metadata, User can also choose to aggregate the gene level matrix by a database ID such as COG, KEGG, DBCAN and PROKKA. '

clustering <- "Clustering with different methods (knn, hdbscan, cut dendrogram) in order to generate new bacterial groups and metadata. Clustering can be performed on the two principal components of the chosen matrix by knn or hdbscan, and cutting the dendrogram created with the chosen matrix "


corepananalysis <- "Overview of the pancore-gene clusters or pancore-BGC clusters by metadata or new metadata generated in clustering section. Pan-genome subsection shows the COG hierarchy functional annotation of the gene clusters present in each genome. User can group genomes by metadata in order to see the variance in functionality and number of genes. Core-genome section presents the functions of the core-gene clusters of the complete dataset as well as the core-gene clusters of metadata groups and their overlap"

statisticalanalysis <- 'Gene level and BGC level statistics between metadata groups, different plots are generated to visualize the statistics results. New metadata groups generated in the clustering section can be used. Statistics can be done at the gene level (CD-HIT , COG, KEGG, DBCAN, PROKKA) and at the BGC level (GCF). Fisher test in combination with basic relative presence statistics are used to study the distribution of the gene clusters among the metadata'

geneclusterdistribution <- "Users can check the distribution of gene clusters of interest by introducing the id of the gene cluster. The Shiny app will generate a PCoA indicating the presence/absence of the gene in the different bacteria"

downloadfasta <- "Users can download the sequences in fasta format of a gene cluster or database id of interest of the complete dataset or specific bacteria"













