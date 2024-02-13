
library(shiny)
library(readr)
library(UpSetR)
library(data.table)
library(dplyr)
library(dbscan)
library(shinydashboard)
library(scales)
library(shinythemes)
library(plotly)
library(DT)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(gplots)
library(gtools)
library(RColorBrewer)
library(viridis)
library(rmarkdown)
library(heatmaply)
library(pals)
library(ape) #Bioconductor manager
library(stringdist)
library(BiocGenerics) #Bioconductor manager
library(Biostrings) #Bioconductor manager
library(tidyr)
library(stringr)
library(vegan)
library(writexl)

source('data_app.R')





# Define UI ----
ui <- navbarPage("bacLIFE", theme = shinytheme("flatly"),
                 tabPanel('INTRODUCTION',
                          sidebarLayout(position = 'left',
                                        sidebarPanel(width = 3,p('Authors: Guillermo Guerrero & Victor Carrion'),
                                                     p('Instituto de Hortofruticultura Subtropical y Mediterránea (IHSM)'),
                                                     p('Universidad de Malaga (UMA)'),
                                                     p('Institute of Biology Leiden (IBL)'),
                                                     p('Leiden University (LU)'),
                                                     p('Netherlands Institute of Ecology (NIOO-KNAW) '),
                                                     imageOutput('image3', height = '60px', width= '200px'),
                                                     
                                        ),
                                        
                                        mainPanel(width = 7,
                                                  h1(tags$b('Welcome to bacLIFE!')),
                                                  p(introduction ),
                                                  imageOutput('image1',  height = '400px', width= '850px'),
                                                  
                                                  h1(tags$b('bacLIFE App')),
                                                  p( style="text-align: justify", bacLIFE_app),
                                                  h3(tags$b('Exploration')),
                                                  p(style="text-align: justify",exploration),
                                                  h3(tags$b('Clustering')),
                                                  p(style="text-align: justify", clustering),
                                                  h3(tags$b('Pancore analysis')),
                                                  p(style="text-align: justify",corepananalysis), 
                                                  h3(tags$b('Statistics')),
                                                  p(style="text-align: justify", statisticalanalysis),
                                                  h3(tags$b('Gene cluster distribution')),
                                                  p(style="text-align: justify",geneclusterdistribution),
                                                  h3(tags$b('Download fasta')),
                                                  p(style="text-align: justify",downloadfasta),
                                                  
                                                  
                                                  p(style="text-align: justify","bacLIFE was created by the CarrionLab of the Institute of Biology Leiden's  (IBL) plant-microbiome interaction department."),
                                                  headerPanel(""),
                                                  headerPanel(""),
                                                  fluidRow(
                                                    column(12,
                                                           div(style="display: inline-block; width: 60%;",imageOutput('image2.1',  height = '100px', width= '600px')),
                                                           div(style="display: inline-block; width: 20%;",imageOutput('image2',  height = '100px', width= '100px')),
                                                           div(style="display: inline-block; width: 20%;",imageOutput('image2.2',  height = '100px', width= '100px'))
                                                    )),
                                                  #column(1,imageOutput('image2.1',  height = '100px', width= '500px')),
                                                  #column(2,imageOutput('image2',  height = '100px', width= '100px')),
                                                  #column(3,imageOutput('image2.2',  height = '100px', width= '100px'))),
                                                  headerPanel(""),
                                                  headerPanel(""),
                                                  br(),
                                                  
                                        ))
                 ),
                 
                 tabPanel('EXPLORATION',
                          tabsetPanel(
                            tabPanel('Gene clusters',
                                     h1(tags$b("Gene clusters")),
                                     sidebarLayout(position = 'left',
                                                   sidebarPanel(width = 3,'Input parameters',
                                                                selectInput('PCA_database', label = "Database:", choices = list('MCL', 'KEGG', 'COG', 'PFAM', 'PROKKA', 'DBCAN')),
                                                                selectInput("PCA_color", label = "Color:", choices = mapping_options_PCA)
                                                   ),
                                                   mainPanel(width = 7,
                                                             tabsetPanel(
                                                               tabPanel('PCoA',
                                                                        tags$label(h3(tags$b('PCoA plot'))),
                                                                        p(style="text-align: justify",'Description: Principal coordinate analysis (PCoA) based on the dice dissimilarity calculated on the absence/presence table of gene clusters. Each point represent one genome and its colored based on the category of the metadata column chosen in the sidebar.'),
                                                                        
                                                                        uiOutput(outputId = "PCAplot",width = '800px',height = '2000px'),
                                                                        tags$label(h3('')),
                                                                        tags$label(h3(tags$b('Permanova test'))),
                                                                        p(style="text-align: justify",'Description: Permanova test to evaluate if there are differences between the groups specified on the category of the metadata column chosen in the sidebar.'),
                                                                        
                                                                        verbatimTextOutput(outputId = "permanova"),
                                                                        br(),
                                                                        br(),
                                                               ),
                                                               tabPanel('Dendrogram',
                                                                        tags$label(h3(tags$b('Dendrogram'))),
                                                                        p(style="text-align: justify",'Description: Dendrogram based on the dice dissimilarity calculated on the absence/presence table of gene clusters. Colors represent the categories of the metadata column chosen in the sidebar.'),
                                                                        
                                                                        plotOutput(outputId = "dendograma",width = '700px',height = '700px'),
                                                                        br(),
                                                                        br(),
                                                               )
                                                             )
                                                   )
                                     )
                                     
                            ),
                            tabPanel('Biosynthetic gene clusters',
                                     h1(tags$b("Biosynthetic gene clusters")),
                                     sidebarLayout(position = 'left',
                                                   sidebarPanel(width = 3,'Input parameters',
                                                                selectInput("PCA_color_BGC", label = "Color:", choices = mapping_options_PCA)
                                                   ),
                                                   mainPanel(width = 7,
                                                             tabsetPanel(
                                                               tabPanel('PCA',
                                                                        p(style="text-align: justify",'Description: Principal Coordinate Analysis based on the dice dissimilarity calculated on the absence/presence table of Biosynthetic Gene Clusters (BGCs). Each point represent one genome and its color the category of the metadata of the metadata column chosen in the sidebar.'),
                                                                        
                                                                        tags$label(h3(tags$b('PCoA plot'))),
                                                                        plotlyOutput(outputId = "PCAplotBGC",width = '800px',height = '500px')),
                                                               br(),
                                                               br(),
                                                               tabPanel('BGCs heatmap',
                                                                        tags$label(h3(tags$b('Heatmap Gene Cluster Families'))),
                                                                        selectInput('bgcheatmap_column', label = 'Side Column color:', choices = mapping_options),
                                                                        
                                                                        
                                                                        p(style="text-align: justify",'Description: Heatmap calculated on the absence/presence table of Biosynthetic Gene Clusters (BGCs). Colors represent the categories of the metadata column chosen in the sidebar.'),
                                                                        
                                                                        plotlyOutput(outputId = "bgcheatmap2",width = '800px',height = '800px'),
                                                                        br(),
                                                                        br(),
                                                                        tags$label(h3(tags$b('Get annotations of a specific Gene Cluster Family (GCF)'))),
                                                                        br(),
                                                                        textInput('GCF_id', 'Select GCF: ', value = ''),
                                                                        actionButton("checkGCF", "Get annotations"),
                                                                        br(),
                                                                        uiOutput(outputId = "GCF_annotations",width = '800px',height = '2000px'),
                                                                        br(),
                                                                        br(),
                                                                        
                                                               ),
                                                               tabPanel('Antismash UI',
                                                                        selectInput('sample4antismash', label = 'Select bacteria to see antismash UI', choices =as.list(unique(mapping_file$Sample) )),
                                                                        actionButton("go_antismashbacteria", "Get antismash"),
                                                                        
                                                                        br(),
                                                                        br(),
                                                                        fluidRow(box(width=NULL, htmlOutput("inc")))
                                                                        
                                                                        
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
                                     h1(tags$b("Gene clusters")),
                                     sidebarLayout(position = 'left',
                                                   sidebarPanel(width = 3,'Input parameters',
                                                                selectInput('clusterclass', label = 'Upset plot data:', choices = list('MCL', 'KEGG', 'COG', 'PFAM', 'PROKKA', 'DBCAN', 'GCF')),
                                                                sliderInput('knumber', label = 'Number of clusters KNN:', min = 1, max = 15, value = 2),
                                                                sliderInput('hdbscannumber', label = 'Minimum points/cluster HDBSCAN:', min = 2, max = 100, value = 2),
                                                                sliderInput('dendonumber', label = 'Number of groups to cut dendrogram:', min = 2, max = 15, value = 2)
                                                                
                                                   ),
                                                   mainPanel(width = 7,
                                                             tabsetPanel(
                                                               tabPanel('KNN clustering',
                                                                        
                                                                        h4(tags$b('Optimal number of clusters')),
                                                                        
                                                                        p(style="text-align: justify",'Description: Within groups sum of squares using 1-15 clusters. This value should be orientative in finding the optimal number of bacterial clusters based in bacLIFE gene cluster absence/presence matrix. '),
                                                                        
                                                                        plotOutput(outputId = "optimalclusters",width = '500px',height = '300px'),
                                                                        h4(tags$b('KNN clusters PCoA plot')),
                                                                        p(style="text-align: justify",'Description: Principal Coordinate Analysis based on the dice dissimilarity calculated on the absence/presence table of gene clusters. Each point represents one genome and its color the cluster membership when K-Nearest Neighbors Algorithm (KNN) is implemented.'),
                                                                        
                                                                        plotlyOutput(outputId = "KNN",width = '800px',height = '500px'),
                                                                        br(),
                                                                        br(),
                                                                        
                                                               ),
                                                               tabPanel('HDBSCAN clustering',
                                                                        h4(tags$b('HDBSCAN clusters PCoA plot')),
                                                                        p(style="text-align: justify",'Description: Principal coordinate analysis based on the dice dissimilarity calculated on the absence/presence table of gene clusters. Each point represent one genome and its color the cluster membership when HDBSCAN clustering is implemented.'),
                                                                        
                                                                        plotlyOutput(outputId = "HDBSCAN",width = '800px',height = '500px'),
                                                                        br(),
                                                                        br(),
                                                                        
                                                               ),
                                                               tabPanel('Dendrogram clustering',
                                                                        h4(tags$b('Dendrogram clusters')),
                                                                        p(style="text-align: justify",'Description: Dendrogram based on the dice dissimilarity calculated on the absence/presence table of gene clusters. Each color represents the cluster membership when the dendrogram is cut.'),
                                                                        
                                                                        plotOutput(outputId = "dendoclust",width = '700px',height = '700px'),
                                                                        br(),
                                                                        br(),
                                                                        
                                                               )
                                                             )))
                                     
                            ),
                            tabPanel("Lifestyle Preservation",
                                     h1(tags$b("Lifestyle preservation")),
                                     sidebarLayout(position = 'left',
                                                   sidebarPanel(width = 3,selectInput('input_column', label = "Select cluster to visualize lifestyle perseverance:", choices = list('clusters_knn', 'clusters_hdbscan', 'clusters_dendrogram')),
                                                                
                                                                
                                                   ),
                                                   mainPanel(width = 7,
                                                             h3(tags$b('Barplots lifestyle preservation per cluster')),
                                                             p(style="text-align: justify",'Description: Barplot showing the proportions of Lifestyles in each of the clusters obtained with any of the three methods seen in the previous tab.'),
                                                             
                                                             plotlyOutput(outputId = "lifestyle_preservation",width = '800px',height = '500px'),
                                                             br(),
                                                             br(),
                                                             
                                                             
                                                   )
                                     ))
                            
                          )
                          
                          
                          
                 ),
                 
                 
                 tabPanel("PAN-CORE ANALYSIS",
                          tabsetPanel(
                            tabPanel('Pan-Genome',  
                                     h1(tags$b("Pan-genome COG profiles")),
                                     sidebarLayout(position = 'left',
                                                   sidebarPanel(width = 3,'Input parameters',
                                                                selectInput('allgenescolumn', label = 'Plot group:', choices = mapping_options),
                                                                #selectInput('allgenesprocess', label = 'Select COG process:', choices = cog_process_options)
                                                   ),
                                                   mainPanel(width = 7,
                                                             
                                                             tags$label(h3(tags$b('COG profiles percentage % '))),
                                                             p(style="text-align: justify",'Description: Barplots showing the COG annotation of all genes of each bacteria. Barplots can be grouped by any metadata column or clusters created in the CLUSTERING section.'),
                                                             
                                                             plotlyOutput(outputId = "allgenescogplot_2",width = '1000px',height = '500px'),
                                                             br(),
                                                             br(),
                                                             
                                                             
                                                             
                                                             
                                                             
                                                   )
                                     )),
                            tabPanel('Core-Genome',  
                                     
                                     h1(tags$b("Core-Genome")),
                                     sidebarLayout(position = 'left',
                                                   sidebarPanel(width = 3,'Input parameters',
                                                                selectInput('upsetclass', label = 'Upset plot data:', choices = list('MCL', 'KEGG', 'COG', 'PFAM', 'PROKKA', 'DBCAN', 'GCF')),
                                                                selectInput('upsetcolumn', label = 'Upset plot group:', choices = mapping_options),
                                                                #selectInput('coreprocess', label = 'Select COG process:', choices = cog_process_options)
                                                   ),
                                                   mainPanel(width = 7,
                                                             tags$label(h3(tags$b('Complete dataset core-genome'))),
                                                             p(style="text-align: justify",'Description: COG functional categories proportions of the core-genome (genes present in >90% of genomes) of the complete dataset.'),
                                                             
                                                             plotlyOutput(outputId = "allcorecog",width = '500px',height = '400px'),
                                                             p(style="text-align: justify",paste0('Number of complete dataset core gene clusters: ', nrow(core_genes))),
                                                             br(),
                                                             p(style="text-align: justify",'Download the core-genome in a table with functional annotations'),
                                                             
                                                             downloadButton("downloadallcore", "Download"),
                                                             br(),
                                                             tags$label(h3(tags$b('Group core genes'))),
                                                             p(style="text-align: justify",'Description: Barplot showing the number of core-genes (genes present in >90% of genomes) of each metadata group.'),
                                                             
                                                             plotlyOutput(outputId = "coreplot"),
                                                             
                                                             uiOutput('table_group4core'),
                                                             br(),
                                                             tags$label(h3(tags$b('Group core COG profiles'))),
                                                             p(style="text-align: justify",'Description: Barplot showing the COG functional categories proportions of each metadata group core-genome.'),
                                                             
                                                             plotlyOutput(outputId = "corecogplot",width = '800px',height = '500px'),
                                                             tags$label(h3(tags$b('Upset plot of core genes overlapping'))),
                                                             p(style="text-align: justify",'Description: The Upset plot displays the overlap of the core-genome among each metadata group.'),
                                                             
                                                             plotOutput(outputId = "Upsetplot",width = '800px',height = '500px'),
                                                             br(),
                                                             br(),
                                                             
                                                             
                                                             
                                                   )
                                     )),
                            
                            tabPanel('Unique Genes',
                                     mainPanel(width = 7,
                                               tags$label(h3(tags$b('Histogram number of unique genes'))),
                                               p(style="text-align: justify",'Description: Histogram displays the frequency of the number of unique genes in the dataset. We consider as a "unique gene" each gene that shows no homology with any other genes of the dataset, hence, being exclusive of one bacterium.'),
                                               
                                               plotlyOutput(outputId = "singletonplot"),
                                               
                                               
                                               tags$label(h3(tags$b('Number of unique genes/sample in this dataset'))),
                                               p(style="text-align: justify",'Description: Table displays the number of unique genes in each bacterium of the dataset.'),
                                               DT::dataTableOutput('table_singletons'),
                                               #box(title = "Number of unique genes/sample in this dataset:", width = NULL, height = 1,solidHeader = T,  status = "primary",div(style = 'height:400px;overflow-y: scroll', tableOutput('table_singletons'))),
                                               tags$label(h3(tags$b('Get unique genes of a specific bacteria'))),
                                               p(style="text-align: justify",'Description: Choose a bacterium to display the unique genes found in that genome and its functional annotations.'),
                                               selectInput('bacteria4uniquegenes', label = "Select bacteria:", choices = as.list(c('All', unique(mapping_file$Sample) ))),
                                               actionButton("go_uniquegenes", "Get unique genes"),
                                               br(),
                                               br(),
                                               DT::dataTableOutput(outputId = "uniquegenestable"),
                                               br(),
                                               #downloadButton("downloadUniques", "Download"),
                                               
                                     )),
                            
                            
                            tabPanel('About',
                                     tags$label(h3(tags$b('COG Categories Abbreviation'))),
                                     DT::dataTableOutput('cogtable'),
                                     
                            )
                            
                            
                          )),
                 
                 
                 
                 tabPanel("STATISTICS",
                          h1(tags$b("Statistics")),
                          sidebarLayout(position = 'left',
                                        sidebarPanel(width = 3,
                                                     selectInput('data_for_stats', label = "Select data to make statistics:", choices = list('MCL', 'KEGG', 'COG', 'PFAM', 'PROKKA', 'DBCAN', 'GCF')),
                                                     selectInput('column_for_stats', label = "Select column to make statistics:", choices = mapping_options),
                                                     uiOutput('inputs1'),
                                                     uiOutput('inputs2'),
                                                     actionButton("go", "Go"),
                                        ),
                                        mainPanel(width = 7,
                                                  tags$label(h3(tags$b('Lifestyle Associated Genes'))),
                                                  p(''),
                                                  p(style="text-align: justify",'Statistics allow users to compare the distributions of genes and BGCs among different metadata groups given by the users or obtained in the CLUSTERING module. Fisher tests together with basic group prevalence statistics allows to find genes that exhibit a distinct pattern of presence for a specific group, while being largely absent in other groups. Please select two groups of interest for the statistical tests, press the button “Go” and proceed to the “Statistics” subsection.'),
                                                  
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
                 tabPanel("GENE DISTRIBUTIONS",
                          sidebarLayout(position = 'left',
                                        sidebarPanel(width = 3,
                                                     h4(tags$b('cluster ID distribution')),
                                                     textInput('cluster2distribution', 'Select cluster: ', value = ''),
                                                     actionButton("checkdistribution", "Get distribution"),
                                                     
                                                     
                                        ),
                                        
                                        mainPanel(width = 7,
                                                  h1(' '),
                                                  h3(tags$b('PCoA with the distribution of the chosen gene cluster')),
                                                  p(style="text-align: justify",'This section allows users to view the distribution of particular gene clusters of interest through a Principal Coordinate Analysis. The cluster identification numbers are denoted as "cluster_xxxxxx" and can be acquired from the statistics section.'),
                                                  
                                                  h1(' '),
                                                  
                                                  br(),
                                                  
                                                  uiOutput(outputId = "PCAplotgenedistribution",width = '800px',height = '2000px'),
                                                  br(),
                                                  br(),
                                        )
                                        
                                        
                                        
                                        
                                        
                          )),
                 tabPanel("DOWNLOAD FASTA",
                          h1(tags$b('Download fasta sequences')),
                          sidebarLayout(position = 'left',
                                        sidebarPanel(width = 3,
                                                     selectInput('data_for_msa', label = "Select bacteria:", choices = as.list(c('All', unique(mapping_file$Sample) ))),
                                                     selectInput('msacolumn', label = 'Select database', choices = list('MCL', 'KEGG', 'COG', 'PFAM', 'PROKKA', 'DBCAN')),
                                                     textInput('msaid', 'Select cluster: ', value = ''),
                                                     #actionButton("gomsa", "Get fasta"),
                                                     downloadButton("downloadFASTA", "Download")
                                                     
                                        ),
                                        
                                        mainPanel(width = 7,
                                                  h1(' '),
                                                  p(style="text-align: justify",'Users can download the sequence of gene clusters/COG/KEGG, PFAM, DBCAN identification numbers of their interest from this section, for all the bacteria in the dataset or specific ones.'),
                                                  h1(' '),
                                                  h3(tags$b('Metadata ')),
                                                  DT::dataTableOutput('map_download'),
                                                  uiOutput('results_fasta_download')
                                                  
                                        )
                          )
                          
                          
                          
                          
                 ),
                 
                 
                 
                 
)







# Define server logic ----
server <- function(input, output){
  
  options(shiny.maxRequestSize=20*1024^2)
  
  
  output$image1 <- renderImage({
    list(src='www/bacLIFE_wokflow.png', height = '400px', width= '850px')
  },  deleteFile = F)
  
  output$image2 <- renderImage({
    list(src='www/Leiden_university.png', height = '100px', width = '100px')
  },  deleteFile = F)
  
  output$image2.1 <- renderImage({
    list(src='www/UMA.png', height = '100px', width = '400px')
  },  deleteFile = F)
  
  output$image2.2 <- renderImage({
    list(src='www/NIOO.png', height = '100px', width = '100px')
  },  deleteFile = F)
  
  
  output$image3 <- renderImage({
    list(src='www/carrionlab.png', height = '60px', width = '200px')
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
    
    fig <- fig %>% layout(title = 'Number of core genes', legend = list(font = list(size = 9), xaxis = list(title = 'Group'), yaxis = list(title = 'Number of Core genes'))
    )
    fig
  })
  
  
  
  output$corecogplot <- renderPlotly({
    column <- input$upsetcolumn
    cog_process <- 'All'
    mapping_file <- new_mapping_file()
    mapping_file <- mapping_file[[1]]
    input <- general_upset_plot(matrix[,1:n_samples], mapping_file, 'full_matrix', column)
    
    plot <- cog_core_groups(matrix, input, cog_extended_annotations, cog_process, n_samples, df_colors)
    ggplotly(plot + ylab('Relative COG abundance %') + xlab('Genome'))
  })
  
  
  df_all_genes_cog <- reactive({
    column <- input$allgenescolumn
    cog_process <- 'All'
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
    ggplotly(all_genes_cog[[2]]+ ylab('Relative COG abundance %') + xlab('Genome'))
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
  
  output$permanova <- renderPrint({
    map = mapping_file
    color_column <- input$PCA_color
    
    reorder_idx <- match(colnames(full_matrix),map$Sample)
    map <- map[reorder_idx,]
    permanova = adonis2(
      formula= distance_matrix ~ map[,color_column], 
      permutations=999,
      method='bray'
    )
    
    print(permanova)
    
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
  
  output$table_singletons <- renderDataTable({
    colnames(singletons_table)[2] = 'Unique genes'
    singletons_table
  })
  
  
  
  
  output$singletonplot <- renderPlotly({
    plot <- ggplot(singletons_table, aes(x=singletons)) + 
      geom_histogram(color="black", fill="white")+ theme_bw() + ylab('Genome count')
    ggplotly(plot)
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
  
  
  vals_uniquegenestable <- eventReactive(input$go_uniquegenes, {
    
    bacteriaid <- input$bacteria4uniquegenes
    
    return(bacteriaid)
  })
  
  
  output$uniquegenestable <-  DT::renderDataTable(server=F,{
    bacteria = vals_uniquegenestable()
    subset_matrix = singletons[,c('clusters', bacteria, 'keggid', 'kegg_description', "cogid", "cog_description", "pfamid",
                                  "pfam_description", "dbcanid", "dbcan_description", "prokkadescription")]
    
    subset_matrix = subset_matrix %>% filter(across(2) > 0)
    substet2print =  DT::datatable(
      subset_matrix,
      caption = paste0("Table with unique genes from: ", bacteria ),
      extensions = c("Buttons"),
      filter = "top",
      options = list(scrollX = TRUE,
                     dom = 'Bfrtip',
                     buttons = list(
                       list(extend = "excel", text = "Download Current Page", filename = "page",
                            exportOptions = list(
                              modifier = list(page = "current")
                            )
                       ), 
                       list(extend = "excel", text = "Download Full Table", filename = "unique_genes",
                            exportOptions = list(
                              modifier = list(page = "all"))))))
    substet2print
    
    
    
    
  })
  
  
  observe({
    if (is.null(input$fastagenomeFile)) return()
    file.copy(input$fastagenomeFile$datapath, "tmp", recursive = T)
    file.rename(paste0('tmp/', tail(strsplit(input$fastagenomeFile$datapath, "/")[[1]], 1)), 'tmp/query_mash.fasta' )
  })
  
  
  
  
  
  input_genome_mash <- reactive({
    req(input$fastagenomeFile)  # Check if a file is selected
    
    cmd = paste0('mash dist mash/reference.msh tmp/query_mash.fasta > tmp/output_mash.txt')
    system(cmd)
    mash_results = read.table('tmp/output_mash.txt')
    mash_results = mash_results[,c(1,3,4,5)]
    colnames(mash_results)= c('bacLIFE genome', 'ANI', 'p-value', 'markers')
    mash_results$ANI = as.numeric(mash_results$ANI)
    mash_results$`p-value` = as.numeric(mash_results$`p-value`)
    mash_results = mash_results[order(mash_results$ANI, decreasing= T),]
    
    return(mash_results)
    
    
  })
  
  output$mash_table <- renderDataTable({
    
    mash_table = input_genome_mash()
    
    
    mash_table$`bacLIFE genome` = str_remove(mash_table$`bacLIFE genome`, './')
    mash_table = merge(names_equivalence, mash_table, by.x = 'bacLIFE_name', by.y = 'bacLIFE genome')
    mash_table$Full_name = str_remove(mash_table$Full_name, '_O.fna')
    mash_table$bacLIFE_name = NULL
    colnames(mash_table)[1] = 'bacLIFE genome'
    
    mash_table <- DT::datatable(
      mash_table,
      extensions = c("Buttons"),
      filter = "top",
      options = list(scrollX = TRUE,
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
    mash_table
    
    
  })
  
  
  output$mashresults <- renderUI({
    req(input$fastagenomeFile)
    tagList(
      h3(tags$b('map2bacLIFE results')),
      br(),
      p(style="text-align: justify",'Average Nucleotide Identity (ANI) similarity of the input genome against the genomes present in this bacLIFE dataset by Mash "https://github.com/marbl/Mash". The lowest the ANI value the more similar are two genomes.'),
      br(),
      DT::dataTableOutput("mash_table"),
      br(),
      br(),
      
    )
    
  })
  
  
  
  
  
  
  observe({
    if (is.null(input$fastaFile)) return()
    file.copy(input$fastaFile$datapath, "tmp", recursive = T)
    file.rename(paste0('tmp/', tail(strsplit(input$fastaFile$datapath, "/")[[1]], 1)), 'tmp/query.fasta' )
  })
  
  
  
  
  
  
  
  input_fasta2 <- reactive({
    req(input$fastaFile)  # Check if a file is selected
    
    file_id = paste0("file_", as.integer(Sys.time()))
    #write.FASTA(input$fastaFile,  paste0( 'tmp/',file_id, '.fasta'))
    cmd = paste0('blastp -query tmp/query.fasta -db blast_db/db -out tmp/', file_id ,'_blast.txt -num_threads 2 -outfmt "6 qseqid sseqid qcovs pident evalue bitscore"')
    system(cmd)
    file_out = paste0('tmp/', file_id ,'_blast.txt')
    descriptions <- matrix[,c('clusters', 'descriptions')]
    descriptions =   separate_rows(descriptions,  descriptions,sep = ",")
    
    results_blast = read.table(file_out, header = F)
    results_blast = results_blast[results_blast$V3 > 50,]
    results_blast = results_blast[results_blast$V5 < 0.05,]
    results_blast2 = results_blast[1:5,]
    results_blast = merge(results_blast, descriptions, by.x = 'V2', by.y = 'descriptions')
    results_blast = results_blast[order(results_blast$V5),]
    results_blast2 = merge(results_blast2, descriptions, by.x = 'V2', by.y = 'descriptions')
    hit = names(which.max(table(results_blast2$clusters)))
    
    return_list = list(results_blast, hit)
    return(return_list)
    
    
    
  })
  
  
  
  
  
  
  output$tableblasthit <- renderDataTable({
    input = input_fasta2()
    blast_results = input[[1]]
    colnames(blast_results) = c('bacLIFE gene', 'query gene', 'Coverage (%)', 'Sequence similarity (%)', 'e-value', 'bitscore', 'bacLIFE gene cluster')
    blast_results <- DT::datatable(
      blast_results,
      extensions = c("Buttons"),
      filter = "top",
      options = list(scrollX = TRUE,
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
    blast_results
    
  })
  
  output$PCAplotgenedistribution2 <- renderUI({
    
    input = input_fasta2()
    hit = input[[2]]
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
    
    
    tagList(
      h3(tags$b('Blast2bacLIFE results')),
      br(),
      p(style="text-align: justify",'Blast2bacLIFE finds the gene cluster with the closest homology to the query gene and performs Principal Coordinate Analysis (PCoA) showing the distribution of the query sequence and its presence in the different genomes included in this dataset.'),
      br(),
      DT::dataTableOutput('tableblasthit'),
      br(),
      p,
      br(),
      
    )
  })
  
  
  
  output$cluster_hit <- DT::renderDataTable(server=F,{
    hit = input$cluster2distribution
    cluster_hit = matrix[matrix$clusters %in% hit,]
    cluster_hit = cluster_hit[,c('clusters', 'keggid', 'kegg_description', 'cogid', 'cog_description', 'pfamid', 'pfam_description', 'dbcanid', 'dbcan_description', 'prokkadescription')]
    
    cluster_hit <- DT::datatable(
      cluster_hit,
      options = list(scrollX = TRUE,
                     dom = 'Bfrtip'))
    cluster_hit
    
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
    
    tagList(
      h3(tags$b('By cluster ID distribution')),
      br(),
      p(style="text-align: justify",'Principal Coordinate Analysis (PCoA) showing the distribution of the query sequence and its presence in the different genomes included in this dataset.'),
      br(),
      DT::dataTableOutput('cluster_hit'),
      br(),
      p,
      br(),
      
    )
    
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
  
  
  output$downloadUniques <- downloadHandler(
    filename = function(){
      paste( "unique_genes", "tsv", sep = ".")
    },
    content = function(file){
      bacteria <- input$bacteria4uniquegenes 
      file2read = paste0('corepan_analysis/singletons_per_sample/',bacteria, '_singletons.txt')
      file_readed = read_delim(file2read, col_names = T, quote = "\"")
      writexl::write_xlsx(file_readed, file)
    } 
  )
  
  
  
  
  results <- reactive({
    input_list <- randomVals()
    mapping_file <- new_mapping_file()
    mapping_file <- mapping_file[[1]]
    #write.table(mapping_file, "mapping_file_stats.txt", row.names = F)
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
    #write.table(result, 'results.txt', row.names = F)
    return(result)
  })
  
  
  output$markdown <- DT::renderDataTable(server=F,{
    results <- results()
    bacteria = input$resultsordered
    if (bacteria %in% 'All'){
      results$descriptions = NULL
      
    }else{
      namess = names_equivalence
      namess$Full_name = str_remove(namess$Full_name, '_O.fna')
      namess$bacLIFE_name = str_remove(namess$bacLIFE_name, '_O.fna')
      
      subnames = namess[namess$Full_name %in% bacteria,]
      patrn = subnames$bacLIFE_name[1]
      results$genome_position <- sapply(results$descriptions, extract_position, patrn)
      genome_position <- as.numeric(results$genome_position)
      results <- separate_rows(results, genome_position, convert = TRUE)
      results$descriptions = NULL
      results = results[results$genome_position > 0,]
    }
    results = results %>% mutate(across(where(is.numeric), ~ round(., 3)))
    results <- DT::datatable(
      results,
      extensions = c("Buttons"),
      filter = "top",
      options = list(scrollX = TRUE,
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
  
  
  
  output$tablemap <- DT::renderDataTable(server=F,{
    table <- new_mapping_file()
    table <- table[[1]]
    DT::datatable(
      table,
      extensions = c("Buttons"),
      filter = "top",
      options = list(
        dom = 'Bfrtip',
        buttons = list(
          list(extend = "excel", text = "Download Table", filename = "data",
               exportOptions = list(
                 modifier = list(page = "all"))))))
  })
  
  output$cogtable <- DT::renderDataTable({
    DT::datatable(cog_extended_annotations2)
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
    selectInput("Group1", "Select Group 1", options)
    
  })
  output$inputs2 <- renderUI({
    groups <- new_mapping_file()
    groups <- groups[[1]]
    column <- input$column_for_stats #Column of mapping file to make the stats with
    db <- input$data_for_stats #data to make the stats
    options <- extract_options(groups, column)
    selectInput("Group2", "Select Group 2", c(options, 'All'))
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
  
  
  
  vals_GCFdistribution <- eventReactive(input$checkGCF, {
    
    gcfid <- input$GCF_id
    
    return(gcfid)
  })
  
  
  
  output$GCF_annotations <- renderUI({
    id = vals_GCFdistribution()
    tagList(
      br(),
      DT::dataTableOutput('BGC_table')
    )
    
  })
  
  output$BGC_table <- DT::renderDataTable(server=F,{
    gcfid = vals_GCFdistribution()
    gcf_subtable= as.data.frame(bgc_descriptions)
    gcf_subtable = gcf_subtable[gcf_subtable$GCF %in% gcfid,]
    gcf_subtable <- DT::datatable(
      gcf_subtable,
      options = list(scrollX = TRUE,
                     dom = 'Bfrtip'))
    gcf_subtable
    
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
      br(),
      br(),
      tags$label(h3(tags$b('Table results'))),
      p(style="text-align: justify",'Description: Table displays the statistics of the performed contrast, with each row representing a gene cluster or database ID based on the chosen options. The columns show the results of a Fisher test and group relative abundance log-fold changes. Additionally, functional annotations of the gene cluster in various databases are included as columns. The table should be filtered base in user interests and research questions '),
      br(),
      selectInput('resultsordered', label = "Select bacteria:", choices = as.list(c('All', unique(mapping_file$Sample) ))),
      p(style="text-align: justify",'If you chose a specific bacteria, an additional column at the end of the table called genome_position displays the position of the gene in the genome. Consecutive numbers indicate consecutive position in the genome.'),
      br(),
      DT::dataTableOutput('markdown'),
      br(),
      tags$label(h3(tags$b('Absence/presence plot all gene clusters'))),
      p(style="text-align: justify",'Description: The Bubble plot displays the distribution of all gene clusters in the selected groups for the contrast. The plot axes represent the relative abundance of each group in the comparison. The bubble size indicates the number of gene clusters that follow that distribution in both groups. Core genes shared between the two groups are located in the top right corner, while unique genes present in one or a few genomes are in the lower left corner. Gene clusters fully present in one group and fully absent in the other group are located in the top left and bottom right corners, respectively.'),
      
      plotlyOutput('abspresplot',width = '700px',height = '500px'),
      #tags$label(h3('Absence/presence plot significant gene clusters')),
      #plotlyOutput('abspresplot2',width = '700px',height = '500px'),
      #tags$label(h3('Volcano plot')),
      #p('Description: Volcano plot of the contrast chosen with genes with a p-value < 0.01 nad logfold > +-2  '),
      
      #plotlyOutput('volcanoplot',width = '700px',height = '500px'),
      br(),
      tags$label(h3(tags$b('COG profiles'))),
      p(style="text-align: justify",'Descripton: Barplot showing the COG profile of gene clusters that show a higher difference between group means than the value specified in the slide bar below. A minimum mean difference of 0.8 will show genes that for example are present in 100% of one group and in less of 20% of the other group (1-0.2 = 0.8) '),
      sliderInput('min_mean_difference', label = 'Minimum mean difference between groups:', min = 0, max = 1, step = 0.05, value = 0.5),
      plotlyOutput('bubbleplot',width = '1000px',height = '700px'),
      br(),
      br(),
      br(),
    )
    
    
  })
  
  output$bubbleplot <- renderPlotly({
    input_list <- randomVals()
    db <- input_list[[3]]
    
    result <- results()
    min_mean_difference <- input$min_mean_difference
    ploty <- COG_significative_genes(result, cog_extended_annotations, min_mean_difference, db)
    ploty
  })
  
  
  
  
  output$downloadallcore <- downloadHandler(
    filename = function(){
      paste( "core_genes", "xlsx", sep = ".")
    },
    content = function(file){
      table = core_genes
      table = table[,c(1, (n_samples+1):ncol(table))]
      table$descriptions = NULL
      table$completeness = NULL
      write_xlsx(table, file)
    } 
  )
  
  
  
  
  output$table_group4core <- renderUI({
    column = input$upsetcolumn
    
    tagList(
      p(style="text-align: justify",'Download the core-genome of a selected group in a table with functional annotations'),
      selectInput('group4core', label = 'Select group', choices = as.list(unique(mapping_file[,column]))),
      downloadButton("downloadgroupcore", "Download"),
      
    )
    
  })
  
  
  reactive_download <- reactive({
    
    metadatagroup = input$group4core
    column = input$upsetcolumn
    db = input$upsetclass
    groups <- new_mapping_file()
    groups <- groups[[1]]
    
    if (db == 'COG'){
      table <- cog_agg_matrix
    }
    if (db == 'KEGG'){
      table <- kegg_agg_matrix
    }
    if (db == 'PFAM'){
      table <- pfam_agg_matrix
    }
    if (db == 'DBCAN'){
      table <- dbcan_agg_matrix
    }
    if (db == 'MCL'){
      table <- matrix[,1:n_samples]
    }
    ##Databases are not implemented yet
    list_cores = general_upset_plot(matrix[,1:n_samples], groups, 'MCL', column)
    core_interest = list_cores[[metadatagroup]]
    todownload = matrix[matrix$clusters %in% core_interest,]
    print(typeof(todownload))
    print(head(todownload))
    todownload = todownload[,c(1, (n_samples+1):ncol(todownload))]
    todownload$descriptions = NULL
    todownload$completeness = NULL
    return(todownload)
  })
  
  
  output$downloadgroupcore <- downloadHandler(
    filename = function(){
      paste( "core_genes", "csv", sep = ".")
    },
    content = function(file){
      todownload = reactive_download()
      write.table(todownload, file, row.names = F, quote = F, sep = ';')
      
    } 
  )
  
  
  vals_antismash_bacteria <- eventReactive(input$go_antismashbacteria, {
    
    bacteriaid <- input$sample4antismash
    
    return(bacteriaid)
  })
  
  
  output$inc<-renderUI({
    
    bacteria = vals_antismash_bacteria()
    names = names_equivalence
    names$Full_name = str_remove( names$Full_name, '_O.fna')
    names$bacLIFE_name = str_remove( names$bacLIFE_name, '.fna')
    names = data.table(names)
    folder = names[names$Full_name %in% bacteria]$bacLIFE_name
    folder_path = paste0('antismash/', folder, '/')
    addResourcePath(folder , folder_path)
    tags$iframe(seamless="seamless",
                src = paste0(folder, "/index.html"), 
                width=1000, height = 800)
    
    
  })
  
  
  
}
# Run the app ----
shinyApp(ui = ui, server = server)

















