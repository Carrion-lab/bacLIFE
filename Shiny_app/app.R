
source('data_app.R')

# Define UI ----
ui <- navbarPage("MICROLIFE", theme = shinytheme("flatly"),
                 tabPanel('INTRODUCTION',
                          sidebarLayout(position = 'left',
                                        sidebarPanel(p('Authors: Guillermo Guerrero; Adrian Pintado & Victor Carrion'),
                                                     p('Institute of Biology Leiden (IBL)'),
                                                     p('Leiden University')
                                        ),
                                        
                                        mainPanel(
                                          h1('Welcome to Microlife!'),
                                          p(introduction ),
                                          imageOutput('image1',  height = '400px', width= '850px'),
                                          
                                          h1('MicroLife App'),
                                          p( 'MicroLife app consists in five sections which offer a complete analysis based in the conservation of genes and biosynthetic genes among genomes. Specific plots are generated in each of these sections helping the user to visualize the large amount of data'),
                                          h3('- Exploration'),
                                          p(exploration),
                                          h3('- Clustering'),
                                          p( clustering),
                                          h3('- Pancore analysis'),
                                          p(corepananalysis), 
                                          h3('- Statistics'),
                                          p( statisticalanalysis),
                                          h3('- Gene cluster distribution'),
                                          p(geneclusterdistribution),
                                          h3('- Download fasta'),
                                          p(downloadfasta),
                                          
                                          
                                          p("MicroLife was created by the CarrionLab of the Institute of Biology Leiden's  (IBL) plant-microbiome interaction department."),
                                          headerPanel(""),
                                          headerPanel(""),
                                          imageOutput('image2',  height = '100px', width= '100px'),
                                          headerPanel(""),
                                          headerPanel(""),
                                          imageOutput('image3', height = '100px', width= '300px'),
                                          headerPanel("")
                                        ))
                 ),
                 
                 tabPanel('EXPLORATION',
                          tabsetPanel(
                            tabPanel('Gene Clusters',
                                     h1("Gene Clusters"),
                                     sidebarLayout(position = 'left',
                                                   sidebarPanel('Input parameters',
                                                                selectInput('PCA_database', label = "Database:", choices = list('MCL', 'KEGG', 'COG', 'PFAM', 'PROKKA', 'DBCAN')),
                                                                selectInput("PCA_color", label = "Color:", choices = mapping_options_PCA)
                                                   ),
                                                   mainPanel(
                                                     tabsetPanel(
                                                       tabPanel('PCA',
                                                                tags$label(h3('PCA plot')),
                                                                uiOutput(outputId = "PCAplot",width = '800px',height = '2000px')),
                                                       tabPanel('Dendrogram',
                                                                tags$label(h3('Dendrogram')),
                                                                plotOutput(outputId = "dendograma",width = '1000px',height = '1000px')
                                                       )
                                                     )
                                                   )
                                     )
                                     
                            ),
                            tabPanel('Biosynthetic gene clusters',
                                     h1("Biosynthetic gene clusters"),
                                     sidebarLayout(position = 'left',
                                                   sidebarPanel('Input parameters',
                                                                selectInput("PCA_color_BGC", label = "Color:", choices = mapping_options_PCA)
                                                   ),
                                                   mainPanel(
                                                     tabsetPanel(
                                                       tabPanel('PCA',
                                                                tags$label(h3('PCA plot')),
                                                                plotlyOutput(outputId = "PCAplotBGC",width = '800px',height = '500px')),
                                                       tabPanel('BGCs Heatmap',
                                                                h1("Biosynthetic Gene Clusters"),
                                                                selectInput('bgcheatmap_column', label = 'Side Column color:', choices = mapping_options),
                                                                tags$label(h3('Heatmap Gene Cluster Family Classes')),
                                                                plotlyOutput(outputId = "bgcheatmap2",width = '800px',height = '800px')
                                                                
                                                       )
                                                     )
                                                   )
                                     )
                            )
                            
                          )
                          
                 ),
                 
                 
                 
                 tabPanel('CLUSTERING',
                          tabsetPanel(
                            tabPanel('Clustering',
                                     h1("Gene Clusters"),
                                     sidebarLayout(position = 'left',
                                                   sidebarPanel('Input parameters',
                                                                selectInput('clusterclass', label = 'Upset plot data:', choices = list('MCL', 'KEGG', 'COG', 'PFAM', 'PROKKA', 'DBCAN', 'GCF')),
                                                                sliderInput('knumber', label = 'Number of clusters KNN:', min = 1, max = 15, value = 2),
                                                                sliderInput('hdbscannumber', label = 'Minimum points/cluster HDBSCAN:', min = 2, max = 100, value = 2),
                                                                sliderInput('dendonumber', label = 'Number of groups to cut dendrogram:', min = 2, max = 15, value = 2)
                                                                
                                                   ),
                                                   mainPanel(
                                                     tabsetPanel(
                                                       tabPanel('KNN clustering',
                                                                
                                                                h4('Optimal number of clusters'),
                                                                plotOutput(outputId = "optimalclusters",width = '500px',height = '300px'),
                                                                h4('KNN clusters plot'),
                                                                plotlyOutput(outputId = "KNN",width = '800px',height = '500px'),
                                                                
                                                       ),
                                                       tabPanel('HDBSCAN clustering',
                                                                h4('HDBSCAN clusters plot'),
                                                                plotlyOutput(outputId = "HDBSCAN",width = '800px',height = '500px')
                                                                
                                                       ),
                                                       tabPanel('Dendrogram clustering',
                                                                h4('Dendrogram clusters'),
                                                                plotOutput(outputId = "dendoclust",width = '1000px',height = '1000px')
                                                                
                                                       )
                                                     )))
                                     
                            ),
                            tabPanel("LIFESTYLE PRESERVATION",
                                     h1("Lifestyle preservation"),
                                     sidebarLayout(position = 'left',
                                                   sidebarPanel(selectInput('input_column', label = "Choose cluster to visualize lifestyle perseverance:", choices = list('clusters_knn', 'clusters_hdbscan', 'clusters_dendrogram')),
                                                                
                                                                
                                                   ),
                                                   mainPanel(
                                                     h3('Barplots lifestyle preservation per cluster'),
                                                     plotlyOutput(outputId = "lifestyle_preservation",width = '800px',height = '500px')
                                                     
                                                     
                                                   )
                                     ))
                            
                          )
                          
                          
                          
                 ),
                 
                 
                 tabPanel("CORE-PAN ANALYSIS",
                          tabsetPanel(
                            tabPanel('PAN-GENOME',  
                                     h1("Pan-genome COG profiles"),
                                     sidebarLayout(position = 'left',
                                                   sidebarPanel('Input parameters',
                                                                selectInput('allgenescolumn', label = 'Plot group:', choices = mapping_options),
                                                                selectInput('allgenesprocess', label = 'Choose COG process:', choices = cog_process_options)
                                                   ),
                                                   mainPanel(
                                                     
                                                     tags$label(h3('COG profiles percent % ')),
                                                     plotlyOutput(outputId = "allgenescogplot_2",width = '1000px',height = '500px')
                                                     
                                                     
                                                     
                                                     
                                                   )
                                     )),
                            tabPanel('CORE-GENOME',  
                                     
                                     h1("Core-genome"),
                                     sidebarLayout(position = 'left',
                                                   sidebarPanel('Input parameters',
                                                                selectInput('upsetclass', label = 'Upset plot data:', choices = list('MCL', 'KEGG', 'COG', 'PFAM', 'PROKKA', 'DBCAN', 'GCF')),
                                                                selectInput('upsetcolumn', label = 'Upset plot group:', choices = mapping_options),
                                                                selectInput('coreprocess', label = 'Choose COG process:', choices = cog_process_options)
                                                   ),
                                                   mainPanel(
                                                     tags$label(h3('Complete dataset core-genome')),
                                                     plotlyOutput(outputId = "allcorecog",width = '500px',height = '400px'),
                                                     tags$label(h3('Number of core genes')),
                                                     plotlyOutput(outputId = "coreplot"),
                                                     tags$label(h3('Core COG profiles')),
                                                     plotlyOutput(outputId = "corecogplot",width = '800px',height = '500px'),
                                                     tags$label(h3('Upset plot of core genes overlapping')),
                                                     plotOutput(outputId = "Upsetplot",width = '800px',height = '500px'),
                                                     
                                                     
                                                     
                                                   )
                                     )),
                            
                            tabPanel('Unique genes',
                                     tags$label(h3('Histogram number of unique genes')),
                                     plotOutput(outputId = "singletonplot"),
                                     box(title = "Number of unique genes/sample in this dataset:", width = NULL, height = 1,solidHeader = T,  status = "primary",div(style = 'height:400px;overflow-y: scroll', tableOutput('table_singletons'))),
                            )
                            
                          )),
                 
                 
                 
                 tabPanel("STATISTICS",
                          h1("Statistics"),
                          sidebarLayout(position = 'left',
                                        sidebarPanel(
                                          selectInput('data_for_stats', label = "Choose data to make statistics:", choices = list('MCL', 'KEGG', 'COG', 'PFAM', 'PROKKA', 'DBCAN', 'GCF')),
                                          selectInput('column_for_stats', label = "Choose column to make statistics:", choices = mapping_options),
                                          uiOutput('inputs1'),
                                          uiOutput('inputs2'),
                                          actionButton("go", "Go"),
                                        ),
                                        mainPanel(
                                          tabsetPanel(
                                            tabPanel("Mapping file for stats:",
                                                     DT::dataTableOutput('tablemap'),
                                                     #box(title = "Mapping file for stats:", width = NULL, height = 1,solidHeader = T,  status = "primary",div(style = 'height:400px;overflow-y: scroll', tableOutput('table')))
                                            ),
                                            tabPanel("Statistics",
                                                     uiOutput('statsresults')
                                                     
                                            )
                                          )
                                          
                                        )
                          )),
                 tabPanel("GENE CLUSTER DISTRIBUTION",
                          h1('Check distribution of a gene cluster'),
                          sidebarLayout(position = 'left',
                                        sidebarPanel(
                                          textInput('cluster2distribution', 'Choose cluster: ', value = ''),
                                          actionButton("checkdistribution", "Obtain distribution"),
                                        
                                          
                                        ),
                                        
                                        mainPanel(
                                          h1(' '),
                                          h3('PCoA with the distribution of the chosen gene cluster'),
                                          uiOutput(outputId = "PCAplotgenedistribution",width = '800px',height = '2000px')
                                          
                                        )
                          )
                          
                          
                          
                          
                 ),
                 tabPanel("DOWNLOAD FASTA",
                          h1('Download fasta sequences'),
                          sidebarLayout(position = 'left',
                                        sidebarPanel(
                                          selectInput('data_for_msa', label = "Choose bacteria:", choices = as.list(c('All', unique(mapping_file$Sample) ))),
                                          selectInput('msacolumn', label = 'Choose database', choices = list('MCL', 'KEGG', 'COG', 'PFAM', 'PROKKA', 'DBCAN')),
                                          textInput('msaid', 'Choose cluster: ', value = ''),
                                          #actionButton("gomsa", "Get fasta"),
                                          downloadButton("downloadFASTA", "Download")
                                          
                                        ),
                                        
                                        mainPanel(
                                          h1(' '),
                                          h1('Metadata '),
                                          DT::dataTableOutput('map_download'),
                                          uiOutput('results_fasta_download')
                                          
                                        )
                          )
                          
                          
                          
                          
                 )
)







# Define server logic ----
server <- function(input, output){
  
  
  
  output$image1 <- renderImage({
    list(src='www/Microlife_wokflow.png', height = '400px', width= '850px')
  },  deleteFile = F)
  
  output$image2 <- renderImage({
    list(src='www/Leiden_university.png', height = '100px', width = '100px')
  },  deleteFile = F)
  
  output$image3 <- renderImage({
    list(src='www/carrionlab.png', height = '100px', width = '300px')
  },  deleteFile = F)
  
  
  df_upset <- reactive({
    mapping_file <- new_mapping_file()
    mapping_file <- mapping_file[[1]]
    group <- input$upsetcolumn
    data <- input$upsetclass
    if (data == 'MCL'){
      plot <- general_upset_plot(matrix[,1:n_samples], mapping_file, 'full_matrix', group)
    }
    if (data == 'KEGG'){
      plot <- general_upset_plot(kegg_agg_matrix, mapping_file, 'kegg', group)
    }
    if (data == 'COG'){
      plot <- general_upset_plot(cog_agg_matrix, mapping_file, 'cog', group)
    }
    if (data == 'PFAM'){
      plot <- general_upset_plot(pfam_agg_matrix, mapping_file, 'pfam', group)
    }
    if (data == 'PROKKA'){
      plot <- general_upset_plot(prokka_agg_matrix, mapping_file, 'prokka', group)
    }
    if (data == 'DBCAN'){
      plot <- general_upset_plot(dbcan_agg_matrix, mapping_file, 'dbcan', group)
    }
    if (data == 'GCF'){
      plot <- general_upset_plot(GCF_matrix, mapping_file, 'GCF', group)
    }
    return(plot)
  }) 
  
  list_upset <- reactive({
    list <- df_upset()
    return(list)
  })
  
  output$coreplot <- renderPlotly({
    input <- list_upset()
    len_vector <- c()
    names_vector <- c()
    for (i in names(input)){
      len <- length(input[[i]])
      len_vector <- c(len_vector, len)
      names_vector <- c(names_vector, i)
    }
    to_plot <- data.frame(Groups = names_vector, Core_genes = len_vector)
    fig <- plot_ly(
      to_plot, x = ~Groups, y = ~Core_genes,
      type = "bar",width = '600px',height = '400px', color = I("black")
    )
    
    fig <- fig %>% layout(title = 'Number of core genes', legend = list(font = list(size = 9))
    )
    fig
  })
  
  
  
  output$corecogplot <- renderPlotly({
    column <- input$upsetcolumn
    cog_process <- input$coreprocess
    mapping_file <- new_mapping_file()
    mapping_file <- mapping_file[[1]]
    input <- general_upset_plot(matrix[,1:n_samples], mapping_file, 'full_matrix', column)
    
    plot <- cog_core_groups(matrix, input, cog_extended_annotations, cog_process, n_samples, df_colors)
    ggplotly(plot)
  })
  
  
  df_all_genes_cog <- reactive({
    column <- input$allgenescolumn
    cog_process <- input$allgenesprocess
    mapping_file <- new_mapping_file()
    mapping_file <- mapping_file[[1]]
    all_genes_cog <- cog_groups_percent(cog_agg_matrix, cog_extended_annotations, mapping_file, column ,cog_process, df_colors)
    return(all_genes_cog)
  })
  output$allgenescogplot_1 <- renderPlotly({
    column <- input$allgenescolumn
    mapping_file <- new_mapping_file()
    mapping_file <- mapping_file[[1]]
    all_genes_cog <- cog_groups(cog_agg_matrix, cog_extended_annotations, mapping_file, column)
    ggplotly(all_genes_cog)
  })
  
  output$allgenescogplot_2 <- renderPlotly({
    all_genes_cog <- df_all_genes_cog()
    ggplotly(all_genes_cog[[2]])
  })
  
  
  output$allcorecog <- renderPlotly({
    cog_core_genes <- cog_core_genes[order(cog_core_genes$x, decreasing = T),]
    
    M <- merge(cog_core_genes,df_colors, by.y = 'processes', by.x = 'COG_processes')
    M <- M[order(M$x, decreasing = T),]
    fig <- plot_ly(M, labels = ~COG_processes, values = ~x, type = 'pie',  width = 800, height = 400,marker = list(colors = M$color) )
    fig <- fig %>% layout(title = 'COG hierarchy core-genome', legend = list(font = list(size = 9)),
                          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
    fig
  })
  
  output$Upsetplot <- renderPlot({
    input <- list_upset()
    upset(fromList(input), order.by = "freq", nsets = length(input), text.scale = 2)
    
  })
  
  df_pca <- reactive({
    data <- input$PCA_database
    if (data == 'MCL'){
      plot_pca <- full_pca
      eig <- full_eig
    }
    if (data == 'KEGG'){
      plot_pca <-kegg_pca
      eig <- kegg_eig
    }
    if (data == 'COG'){
      plot_pca <-cog_pca
      eig <- cog_eig
    }
    if (data == 'PFAM'){
      plot_pca <-pfam_pca
      eig <- pfam_eig
    }
    if (data == 'PROKKA'){
      plot_pca <-prokka_pca
      eig <- prokka_eig
    }
    if (data == 'DBCAN'){
      plot_pca <-dbcan_pca
      eig <- dbcan_eig
    }
    if (data == 'GCF'){
      plot_pca <-GCF_pca
      eig <- GCF_eig
    }
    return_list <- list(plot_pca,eig)
    return(return_list)
    
  })
  explorative_dendogram <- reactive({
    data <- input$PCA_database
    if (data == 'MCL'){
      dend <- full_dend
    }
    if (data == 'KEGG'){
      dend <- kegg_dend
    }
    if (data == 'COG'){
      dend <- cog_dend
    }
    if (data == 'PFAM'){
      dend <- pfam_dend
    }
    if (data == 'PROKKA'){
      dend <- prokka_dend
    }
    if (data == 'DBCAN'){
      dend <- dbcan_dend
    }
    if (data == 'GCF'){
      dend <- GCF_dend
    }
    return(dend)
  })
  
  ranges2 <- reactiveValues(x = NULL, y = NULL)
  output$PCAplotBGC <- renderPlotly({
    plot_pca <-GCF_pca
    eig <- GCF_eig
    input_list <- list(plot_pca, eig )
    color_column <- input$PCA_color_BGC
    plot <- ggplot(input_list[[1]], aes_string(x = 'X', y = 'Y', color = color_column, label= 'Sample'))  + xlab(paste0('PC1 - ', round(input_list[[2]][1],2), '%', sep='')) + 
      ylab(paste0('PC2 - ', round(input_list[[2]][2],2), '%', sep='')) + geom_point(aes_string(color = color_column)) + theme_bw() +  scale_color_manual(values = as.character(glasbey(n = nrow(input_list[[1]]))))
    ggplotly(plot)
    
  })
  
  output$PCAplot <- renderUI({
    input_list <- df_pca()
    
    color_column <- input$PCA_color
    
    m <- highlight_key(input_list[[1]])
    
    plot <- ggplot(m, aes_string(x = 'X', y = 'Y', color = color_column, label= 'Sample'))  + xlab(paste0('PC1 - ', round(input_list[[2]][1],2), '%', sep='')) + 
      ylab(paste0('PC2 - ', round(input_list[[2]][2],2), '%', sep='')) + geom_point(aes_string(color = color_column)) + theme_bw() +  scale_color_manual(values = as.character(glasbey(n = nrow(input_list[[1]]))))
    gg <- highlight(ggplotly(plot, height = 500), 'plotly_selected')
    p <- crosstalk::bscols(gg, DT::datatable(m), widths = 10,  device = 'lg')
    p
  })
  
  output$dendograma <- renderPlot({
    hc <- explorative_dendogram()
    dend <- as.phylo(hc)
    groups <- mapping_file
    rownames(groups) <- groups$Sample
    groups <- groups[dend$tip.label,]
    groups <- data.table(groups)
    cols <- assign_colors(groups[[input$PCA_color]]  )
    unique_values <- unique(groups[[input$PCA_color]])
    
    unique_colors <- glasbey(n= length(unique_values))
    co <- rep("black", Nedge(dend))
    for (i in 1:length(unique_values)){
      edges <- which.edge(dend, groups[groups[[input$PCA_color]] %in% unique_values[i]]$Sample)
      co[edges] <- unique_colors[i]
    }
    
    dend$tip.label[] <- ""
    
    plot(dend, type = "fan", no.margin = TRUE, edge.col = co)
    
    
    legend("topright",
           legend = unique_values,
           fill = unique_colors,       # Color of the squares
           border = "black")
    
  })
  
  df_clustering <- reactive({
    data <- input$clusterclass
    if (data == 'MCL'){
      plot_pca <- full_pca
    }
    if (data == 'KEGG'){
      plot_pca <-kegg_pca
    }
    if (data == 'COG'){
      plot_pca <-cog_pca
    }
    if (data == 'PFAM'){
      plot_pca <-pfam_pca
    }
    if (data == 'PROKKA'){
      plot_pca <-prokka_pca
    }
    if (data == 'DBCAN'){
      plot_pca <-dbcan_pca
    }
    if (data == 'GCF'){
      plot_pca <-GCF_pca
    }
    return(plot_pca)
    
  })
  
  db_dendogram<- reactive({
    data <- input$clusterclass
    if (data == 'MCL'){
      dend <- full_dend
    }
    if (data == 'KEGG'){
      dend <- kegg_dend
    }
    if (data == 'COG'){
      dend <- cog_dend
    }
    if (data == 'PFAM'){
      dend <- pfam_dend
    }
    if (data == 'PROKKA'){
      dend <- prokka_dend
    }
    if (data == 'DBCAN'){
      dend <- dbcan_dend
    }
    if (data == 'GCF'){
      dend <- GCF_dend
    }
    return(dend)
  })
  
  
  output$optimalclusters <- renderPlot({
    pca_data <- df_clustering()
    wss <- (nrow(pca_data[,2:3])-1)*sum(apply(pca_data[,2:3],2,var))
    for (i in 2:15) wss[i] <- sum(kmeans(pca_data[,2:3],
                                         centers=i)$withinss)
    plot(1:15, wss, type="b", xlab="Number of Clusters",
         ylab="Within groups sum of squares")
  })
  
  
  output$KNN <- renderPlotly({
    knn<- new_mapping_file()
    knn<- knn[[2]]
    knn$Sample <- paste0(knn$Sample, ' (', knn$Lifestyle,')')
    plot <- ggplot(knn, aes(x = X, y = Y, color = cluster, label = Sample)) + xlab("PC1")+
      ylab("PC2") + geom_point(aes(color = cluster)) + theme_bw()  +  scale_color_manual(values = as.character(glasbey(n = nrow(knn))))
    ggplotly(plot)
  })
  
  output$HDBSCAN <- renderPlotly({
    hdbscan<- new_mapping_file()
    hdbscan<- hdbscan[[3]]
    hdbscan$Sample <- paste0(hdbscan$Sample, ' (', hdbscan$Lifestyle,')')
    plot <- ggplot(hdbscan, aes(x = X, y = Y, color = cluster, label = Sample)) + xlab("PC1")+
      ylab("PC2") + geom_point(aes(color = cluster)) + theme_bw() +  scale_color_manual(values = as.character(glasbey(n = nrow(hdbscan))))
    ggplotly(plot)
  })
  
  new_mapping_file <- reactive({
    
    knn_return<- clustering_knn(df_clustering() , input$knumber)
    knn <- knn_return
    knn$cluster <- paste0('cluster_', knn$cluster)
    colnames(knn)[grep("cluster", colnames(knn))] <- 'clusters_knn'
    knn<- knn[,c('Sample', 'clusters_knn')]
    
    new_df <- merge(mapping_file,knn , by= 'Sample')
    
    
    hdbscan_return<- clustering_hdbscan(df_clustering() , input$hdbscannumber)
    hdbscan <- hdbscan_return
    hdbscan$cluster <- paste0('cluster_', hdbscan$cluster)
    colnames(hdbscan)[grep("cluster", colnames(hdbscan))] <- 'clusters_hdbscan'
    hdbscan<- hdbscan[,c('Sample', 'clusters_hdbscan')]
    new_df <- merge(new_df, hdbscan, by= 'Sample')
    #new_df <- new_df[new_df$cluster_hdbscan != 'Outliers' ]
    
    dendograma <- db_dendogram()
    dend <- as.phylo(dendograma)
    clusters_dendogram <- data.frame(clusters_dendogram = cutree(dendograma, input$dendonumber))
    clusters_dendogram$clusters_dendogram <- paste0('cluster_', clusters_dendogram$clusters_dendogram)
    clusters_dendogram$Sample <- rownames(clusters_dendogram) 
    new_df <- merge(new_df, clusters_dendogram, by= 'Sample')
    
    return_list <- list(new_df, knn_return, hdbscan_return, dend)
    return(return_list)
  })
  
  output$dendoclust <- renderPlot({
    input<- new_mapping_file()
    groups <- input[[1]]
    dend <- input[[4]]
    rownames(groups) <- groups$Sample
    groups <- groups[dend$tip.label,]
    groups <- data.table(groups)
    cols <- assign_colors(groups[,'clusters_dendogram'])
    
    unique_values <- unique(groups[,get('clusters_dendogram')])
    
    unique_colors <- glasbey(n= length(unique_values))
    co <- rep("black", Nedge(dend))
    for (i in 1:length(unique_values)){
      edges <- which.edge(dend, groups[groups[,get('clusters_dendogram')] %in% unique_values[i]]$Sample)
      co[edges] <- unique_colors[i]
    }
    
    dend$tip.label[] <- ""
    
    plot(dend, type = "fan", no.margin = TRUE, edge.col = co)
    
    
    legend("topright",
           legend = unique_values,
           fill = unique_colors,       # Color of the squares
           border = "black")
    
  })
  
  
  
  output$lifestyle_preservation <- renderPlotly({
    column <- input$input_column
    mapping_file <- new_mapping_file()
    mapping_file <- mapping_file[[1]]
    df <- mapping_file %>% group_by(get(column), Lifestyle) %>% summarise(n = n()) %>% mutate(freq = n / sum(n))
    colnames(df)[1] <- column
    df$freq <- 100*df$freq
    df$freq <- round(df$freq, digits=2)
    #df$freq <- paste0(as.character(df$freq), '%')
    plot <- ggplot(df, aes(x=Lifestyle, y = freq, fill = Lifestyle)) + 
      geom_bar(stat="identity") + 
      facet_wrap(as.formula(paste("~", column))) +  
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(axis.text.x=element_blank()) + 
      geom_text(aes(label= paste0(as.character(freq), '%')))  +  scale_fill_manual(values = as.character(polychrome(n = nrow(df)))) + theme(strip.background = element_rect(colour = 'black', fill = 'white'))
    ggplotly(plot)
  })
  
  
  output$table_singletons <- renderTable({
    singletons_table
  })
  
  output$singletonplot <- renderPlot({
    qplot(singletons_table$singletons, geom="histogram", xlab = "Number of singleton genes", ylab= 'Number of samples' ,main = "Histogram of number of singletons",fill=I("blue"), bins = n_samples/2)
  })
  
  
  output$download <- downloadHandler(
    filename = function(){"mapping_file_stats.txt"}, 
    content = function(fname){
      write.table(new_mapping_file(), fname, row.names = F)
    }
  )
  
  
  randomVals <- eventReactive(input$go, {
    
    groups <- new_mapping_file()
    groups <- groups[[1]]
    column <- input$column_for_stats #Column of mapping file to make the stats with
    db <- input$data_for_stats #data to make the stats
    group1 <- input$Group1
    group2 <- input$Group2
    list<- list(groups,column, db, group1, group2)
    return(list)
  })
  
  vals_gomsa <- eventReactive(input$gomsa, {
    mapping_file <- new_mapping_file()
    mapping_file <- mapping_file[[1]]
    id <- input$msaid
    bacteria <- input$data_for_msa 
    db <- input$msacolumn
    return_list <- list(id, mapping_file, db, bacteria)
    return(return_list)
  })
  
  
  vals_checkdistribution <- eventReactive(input$checkdistribution, {
    
    clusterid <- input$cluster2distribution
    
    return(clusterid)
  })
  
  
  output$PCAplotgenedistribution <- renderUI({
    hit = vals_checkdistribution()
    cluster_hit = matrix[matrix$clusters %in% hit,]
    cluster_hit = cluster_hit[,2:n_samples]
    copy_pcoa_data = full_pca
    new_col = as.data.frame(t(cluster_hit))
    new_col$Sample = colnames(cluster_hit)
    colnames(new_col)[1]= 'Presence' 
    new_col = new_col[,c(2,1)]
    copy_pcoa_data = data.table(merge(full_pca, new_col, by = 'Sample'))
    #copy_pcoa_data[copy_pcoa_data$Presence > 1 ]$Presence  = 1
    copy_pcoa_data$Presence = as.factor(copy_pcoa_data$Presence)
    
    
    
    color_column <- 'Presence'
    
    m <- highlight_key(copy_pcoa_data)
    
    plot <- ggplot(m, aes_string(x = 'X', y = 'Y', color = color_column, label= 'Sample'))  + xlab(paste0('PC1 - ', round(full_eig[1],2), '%', sep='')) + 
      ylab(paste0('PC2 - ', round(full_eig[2],2), '%', sep='')) + geom_point(aes_string(color = color_column)) + theme_bw() +  scale_color_manual(values = as.character(glasbey(n = nrow(copy_pcoa_data))))
    gg <- highlight(ggplotly(plot, height = 500), 'plotly_selected')
    p <- crosstalk::bscols(gg, DT::datatable(m), widths = 10,  device = 'lg')
    p
  })
  
  
  
  
  
  fasta_process <- reactive({
    input_list <- vals_gomsa()
    returned_item = extract_sequences_fasta(matrix, all_proteins, input_list[[3]], input_list[[4]], input_list[[1]], names_equivalence)
    
    return(returned_item)
  })
  
  output$downloadFASTA <- downloadHandler(
    filename = function(){
      paste( "sequence_file", "fasta", sep = ".")
    },
    content = function(file){
      mapping_file <- new_mapping_file()
      mapping_file <- mapping_file[[1]]
      id <- input$msaid
      bacteria <- input$data_for_msa 
      db <- input$msacolumn
      input_list1 <- list(id, mapping_file, db, bacteria)
      
      input_list2 = extract_sequences_fasta(matrix, all_proteins, input_list1[[3]], input_list1[[4]], input_list1[[1]], names_equivalence)
      ape::write.FASTA(input_list2[[2]], file)
    } 
  )
  
  
  
  
  results <- reactive({
    input_list <- randomVals()
    mapping_file <- new_mapping_file()
    mapping_file <- mapping_file[[1]]
    write.table(mapping_file, "mapping_file_stats.txt", row.names = F)
    if (input_list[[3]] == 'GCF'){
      withProgress(message = 'Rendering, please wait!', {
        result <- stats_gcf(input_list, GCF_matrix)
      })
    }
    else if (input_list[[3]] == 'MCL'){
      withProgress(message = 'Rendering, please wait!', {
        result <- stats(input_list, clean_matrix)
      })
      
    }else{
      withProgress(message = 'Rendering, please wait!', {
        result <- stats(input_list, matrix)
      }) }
    write.table(result, 'results.txt', row.names = F)
    return(result)
  })
  
  
  output$markdown <- DT::renderDataTable({
    results <- results()
    results$descriptions = NULL
    results <- DT::datatable(
      results,
      extensions = c("Buttons"),
      options = list(
        dom = 'Bfrtip',
        buttons = list(
          list(extend = "excel", text = "Download Current Page", filename = "page",
               exportOptions = list(
                 modifier = list(page = "current")
               )
          ),
          list(extend = "excel", text = "Download Full Results", filename = "data",
               exportOptions = list(
                 modifier = list(page = "all"))))))
    results
  })
  
  
  
  output$tablemap <- DT::renderDataTable({
    table <- new_mapping_file()
    table <- table[[1]]
    DT::datatable(
      table,
      extensions = c("Buttons"),
      options = list(
        dom = 'Bfrtip',
        buttons = list(
          list(extend = "excel", text = "Download Table", filename = "data",
               exportOptions = list(
                 modifier = list(page = "all"))))))
  })
  
  output$map_download <- DT::renderDataTable({
    mapping_file <- new_mapping_file()
    mapping_file <- mapping_file[[1]]
    mapping_file
  })
  
  output$inputs1 <- renderUI({
    groups <- new_mapping_file()
    groups <- groups[[1]]
    
    column <- input$column_for_stats #Column of mapping file to make the stats with
    db <- input$data_for_stats #data to make the stats
    options <- extract_options(groups, column)
    selectInput("Group1", "Choose Group 1", options)
    
  })
  output$inputs2 <- renderUI({
    groups <- new_mapping_file()
    groups <- groups[[1]]
    column <- input$column_for_stats #Column of mapping file to make the stats with
    db <- input$data_for_stats #data to make the stats
    options <- extract_options(groups, column)
    selectInput("Group2", "Choose Group 2", c(options, 'All'))
  })
  
  output$abspresplot <- renderPlotly({
    input_list <- randomVals()
    groups <- input_list[[1]]
    column <- input_list[[2]] #Column of mapping file to make the stats with
    group1 <- input_list[[4]]
    group2 <- input_list[[5]]
    plot <- super_abs_pres_plot(matrix[,c(1:n_samples)], groups, group1, group2, column)
    plot 
  })
  
  
  output$abspresplot2 <- renderPlotly({
    input_list <- randomVals()
    db <- input_list[[3]]
    result<- results()
    plot <- abs_pres_plot(result, db)
    ggplotly(plot) 
  })
  
  output$volcanoplot <- renderPlotly({
    input_list <- randomVals()
    db <- input_list[[3]]
    result<- results()
    plot <- vulcano_plot(result, db)
    ggplotly(plot)
  })
  
  output$bgcheatmap <- renderPlotly({
    groups <- new_mapping_file()
    groups <- groups[[1]]
    column <- input$bgcheatmap_column
    
    heatmap_bgc <- heatmap_gcf_function(GCF_matrix, groups, column)
    heatmap_bgc 
  })
  
  output$bgcheatmap2 <- renderPlotly({
    groups <- new_mapping_file()
    groups <- groups[[1]]
    column <- input$bgcheatmap_column
    
    heatmap_bgc <- heatmap_gcf_function2(GCF_matrix, groups, column)
    heatmap_bgc 
  })
  
  output$results_fasta_download <- renderUI({
    input_list <- vals_gomsa()
    file_download = fasta_process()
    tagList(
      tags$label(h3(paste0('Fasta files saved as ',file_download[[1]] ))),
    )
    
  })
  
  output$statsresults <- renderUI({
    input_list <- randomVals()
    tagList(
      tags$label(h3('Table results')),
      DT::dataTableOutput('markdown'),
      tags$label(h3('Absence/presence plot all gene clusters')),
      plotlyOutput('abspresplot',width = '700px',height = '500px'),
      tags$label(h3('Absence/presence plot significant gene clusters')),
      plotlyOutput('abspresplot2',width = '700px',height = '500px'),
      tags$label(h3('Volcano plot')),
      plotlyOutput('volcanoplot',width = '700px',height = '500px'),
      tags$label(h3('COG profiles')),
      sliderInput('min_mean_difference', label = 'Minimum mean difference between groups:', min = 0, max = 1, step = 0.05, value = 0.5),
      plotlyOutput('bubbleplot',width = '1000px',height = '700px'))
    
  })
  
  output$bubbleplot <- renderPlotly({
    input_list <- randomVals()
    db <- input_list[[3]]
    
    result <- results()
    min_mean_difference <- input$min_mean_difference
    ploty <- COG_significative_genes(result, cog_extended_annotations, min_mean_difference, db)
    ploty
  })
  
}
# Run the app ----
shinyApp(ui = ui, server = server)



















