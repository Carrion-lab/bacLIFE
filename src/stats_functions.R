####################################
###     Guilermo Guerrero        ###
### G.GuerreroEgido@nioo.knaw.nl ###
###           2021               ###
###     Leiden University        ###
###                              ###
####################################




library(dplyr)
library(data.table)
library(viridis)
library(ggplot2)
library(gtools)
library(rmarkdown)
library(reshape2)
#source('corepan_functions.R')

#Remove rows with all colums same value 
#Remove rows with all colums same value 
remove_constant_columns <- function(matrix){
  n_samples <- ncol(matrix)
  keep <- apply(matrix[2:n_samples], 1, function(x) length(unique(x[!is.na(x)])) != 1)
  new_matrix <-  matrix[keep,]
  return(new_matrix)
}


###FISHER TEST FUNCTION###
fisher_test <- function(matrix, groups){
  #prepare data for fisher test
  matrix<- t(matrix)
  colnames(matrix)<- matrix[1,]
  matrix <- matrix[2:length(matrix[,1]),]
  matrix <- as.data.frame(matrix)
  matrix$Sample <- rownames(matrix)
  matrix <- merge(matrix, groups, by= 'Sample')
  
  
  n_variables <- grep("Lifestyle", colnames(matrix)) - 1  
  matrix[2:n_variables] <- apply(matrix[2:n_variables], 2, as.numeric) #convert to numeric the values
  matrix$Lifestyle <- as.factor(matrix$Lifestyle)
  matrix[2:n_variables] <- apply(matrix[2:n_variables], 2, function(x) ifelse(x > 1, 1, x)) #Convert numbers higher than 1 to 1 (categorical data)
  matrix <- Filter(function(x) length(unique(x))> 1, matrix) #Remove constant columns
  n_variables <- grep("Lifestyle", colnames(matrix)) - 1  
  matrix$Lifestyle <- as.factor(unclass(as.factor(matrix$Lifestyle)))
  
  ##Fisher test
  pvalues <- c()
  estimates <- c()
  for (i in 2:n_variables){
    col <- matrix[,i]
    tbl <- table(matrix$Lifestyle, col)
    ftest<- fisher.test(tbl)
    e <- ftest$estimate
    p <- ftest$p.value
    pvalues <- c(pvalues, p)
    estimates <- c(estimates, e)
  }
  
  fisher_pvalues <- data.frame(id = colnames(matrix)[2:n_variables], pvalue_fisher = pvalues)
  return(fisher_pvalues)
}



###LOGISTIC REGRESSION FUNCTION###
logistic_regression <- function(matrix, groups){
  #prepare data for logistic regression
  matrix<- t(matrix)
  colnames(matrix)<- matrix[1,]
  matrix <- matrix[2:length(matrix[,1]),]
  matrix <- as.data.frame(matrix)
  matrix$Sample <- rownames(matrix)
  matrix <- merge(matrix, groups, by= 'Sample')
  group <- unique(matrix$Lifestyle)
  group1 <- group[1]
  group2 <- group[2]
  
  
  n_variables <- grep("Lifestyle", colnames(matrix)) - 1  
  matrix[2:n_variables] <- apply(matrix[2:n_variables], 2, as.numeric)
  matrix$Lifestyle <- as.factor(matrix$Lifestyle)
  matrix$Lifestyle <- as.factor(unclass(as.factor(matrix$Lifestyle)))
  
  
  
  #Logistic regression
  pvalues <- c()
  foldchanges <- c()
  lifestyles <- unique(as.character(matrix$Lifestyle))
  for (i in 2:n_variables){
    sub <- matrix[,c(i,n_variables+1)]
    model <- glm( Lifestyle ~ matrix[,i], data= matrix , family = binomial)
    coef <- summary(model)$coef
    p <- coef[8]
    pvalues <- c(pvalues, p)
    
    sub<- data.table(sub)
    colnames(sub)[1]<-'id'
    
    group1_mean <- mean(sub[sub$Lifestyle ==lifestyles[1],]$id)
    group2_mean <- mean(sub[sub$Lifestyle ==lifestyles[2],]$id)
    
    foldchange <- foldchange(group1_mean,group2_mean)
    foldchanges <- c(foldchanges, foldchange)
    
    
  }
  log2fold <- foldchange2logratio(foldchanges, base=2)
  logit_coef <- data.frame(id = colnames(matrix)[2:n_variables], pvalue_log = pvalues, log2fold = log2fold)
  return(logit_coef)
}


###Calculate fold change (need to be modified so if one value is 0 change it to 0.01)
fold_change <- function(matrix, groups){
  #prepare data for logistic regression
  matrix<- t(matrix)
  colnames(matrix)<- matrix[1,]
  matrix <- matrix[2:length(matrix[,1]),]
  matrix <- as.data.frame(matrix)
  matrix$Sample <- rownames(matrix)
  matrix <- merge(matrix, groups, by= 'Sample')
  group <- unique(matrix$Lifestyle)
  group1 <- group[1]
  group2 <- group[2]
  
  
  n_variables <- grep("Lifestyle", colnames(matrix)) - 1  
  matrix[2:n_variables] <- apply(matrix[2:n_variables], 2, as.numeric)
  matrix$Lifestyle <- as.factor(matrix$Lifestyle)
  matrix$Lifestyle <- as.factor(unclass(as.factor(matrix$Lifestyle)))
  
  
  foldchanges <- c()
  lifestyles <- unique(as.character(matrix$Lifestyle))
  for (i in 2:n_variables){
    sub <- matrix[,c(i,n_variables+1)]
    
    sub<- data.table(sub)
    colnames(sub)[1]<-'id'
    
    group1_mean <- mean(sub[sub$Lifestyle ==lifestyles[1],]$id)
    group2_mean <- mean(sub[sub$Lifestyle ==lifestyles[2],]$id)
    
    foldchange <- foldchange(group1_mean,group2_mean)
    foldchanges <- c(foldchanges, foldchange)
    
    
  }
  log2fold <- foldchange2logratio(foldchanges, base=2)
  log2 <- data.frame(id = colnames(matrix)[2:n_variables], log2fold = log2fold)
  return(log2)
}


##Extract relationships between groups of each gene/id/...
extract_relationships <- function(matrix, groups){
  matrix<- t(matrix)
  colnames(matrix)<- matrix[1,]
  matrix <- matrix[2:length(matrix[,1]),]
  matrix <- as.data.frame(matrix)
  matrix$Sample <- rownames(matrix)
  matrix <- merge(matrix, groups, by= 'Sample')
  
  
  n_variables <- grep("Lifestyle", colnames(matrix)) - 1  
  matrix[2:n_variables] <- apply(matrix[2:n_variables], 2, as.numeric)
  matrix[2:n_variables] <- apply(matrix[2:n_variables], 2, function(x) ifelse(x > 1, 1, x)) #Convert numbers higher than 1 to 2 (categorical data)
  matrix <- Filter(function(x) length(unique(x))> 1, matrix) #Remove constant columns
  n_variables <- grep("Lifestyle", colnames(matrix)) - 1  
  
  lifestyles <- unique(as.character(matrix$Lifestyle))
  
  group1 <- c()
  group2 <- c()
  
  for (i in 2:n_variables){
    col <- matrix[,i]
    tbl <- table(matrix$Lifestyle, col)
    percentages <- percent_presence(tbl, lifestyles)
    group1 <- c(group1, percentages[[1]])
    group2 <- c(group2, percentages[[2]])
    
  }
  
  relationships_df <- data.frame(id = colnames(matrix)[2:n_variables], group1 = group1, group2 = group2)
  colnames(relationships_df)[2] <- lifestyles[1]
  colnames(relationships_df)[3] <- lifestyles[2]
  
  return(relationships_df)
}

###Extract significant genes and calculate qp adjusted valusevalues (pvalues corrected for a FDR of 0.05) for cd-hit stats
extract_significant_genes_cdhit <- function(matrix, db){
  matrix <- matrix[!is.na(matrix$pvalue_fisher),]
  matrix <- data.table(matrix)
  
  matrix[!between(matrix$pvalue_fisher, 0, 1)]$pvalue_fisher <- 1
  
  matrix <- as.data.frame(matrix)
  
  
  matrix$pvalue_adj_fisher <-  p.adjust(matrix$pvalue_fisher, method = "BH", n=length(matrix$pvalue_fisher))
  significant_genes <- matrix#[matrix$qvalue_fisher < 0.05,]
  #significant_genes$pvalue_log <- NULL
  #significant_genes <- significant_genes[,c(1,2,5,3,4)]
  return(significant_genes)
}
  
  
  
###Extract significant genes and calculate qp adjusted valusevalues (pvalues corrected for a FDR of 0.05) for non cd-hit stats
  
extract_significant_genes <- function(matrix, db){
  matrix <- matrix[!is.na(matrix$pvalue_fisher),]
  matrix <- matrix[!is.na(matrix$pvalue_log),]
  matrix <- data.table(matrix)
  
  matrix[matrix$pvalue_log > 0.9,]$pvalue_log <- matrix[matrix$pvalue_log > 0.9,]$pvalue_fisher
  matrix[matrix$log2fold == Inf,]$pvalue_log <- matrix[matrix$log2fold == Inf,]$pvalue_fisher
  matrix[matrix$log2fold == -Inf,]$pvalue_log <- matrix[matrix$log2fold == -Inf,]$pvalue_fisher
  
  matrix[!between(matrix$pvalue_fisher, 0, 1)]$pvalue_fisher <- 1
  matrix[!between(matrix$pvalue_log, 0, 1)]$pvalue_log <- 1
  
  matrix <- as.data.frame(matrix)
  
  
  matrix$pvalue_adj_fisher <-  p.adjust(matrix$pvalue_fisher, method = "bonferroni", n=length(matrix$pvalue_fisher))
  matrix$pvalue_adj_log <-  p.adjust(matrix$pvalue_log, method = "bonferroni", n=length(matrix$pvalue_log))
  significant_genes <- matrix#[matrix$qvalue_fisher < 0.05,]
  #significant_genes$pvalue_log <- NULL
  #significant_genes <- significant_genes[,c(1,2,5,3,4)]
  return(significant_genes)
}


###Extract group gene presence 
percent_presence <- function(tbl, lifestyles){
  group1 <- rownames(tbl)[which(rownames(tbl) == lifestyles[1])]
  percent_group1<- tbl[group1,colnames(tbl)== 1] / (tbl[group1,colnames(tbl)== 0] + tbl[group1,colnames(tbl)== 1])
  
  group2 <- rownames(tbl)[which(rownames(tbl) == lifestyles[2])]
  percent_group2<- tbl[group2,colnames(tbl)== 1] / (tbl[group2,colnames(tbl)== 0] + tbl[group2,colnames(tbl)== 1])
  
  return_list<- list(percent_group1, percent_group2)
  return(return_list)
}








###STATSLOOP cd-hit
superstats_cdhit <- function(data_stats,annotations, column, new_mapping){
    
    full_matrix<- remove_constant_columns(data_stats)
    full_fisher <- fisher_test(data_stats, new_mapping)
    
    full_logfold <- fold_change(data_stats,new_mapping)
    full_relationships <- extract_relationships(data_stats, new_mapping)
    
    
    full_matrix_stats <- merge(full_fisher, full_logfold, by = 'id', all.y = T)
    full_matrix_stats <- merge(full_matrix_stats, full_relationships, by= 'id', all = T)
    full_significant_genes <- extract_significant_genes_cdhit(full_matrix_stats, 'Full_matrix')
    full_significant_genes_annotated <- merge(full_significant_genes, annotations,by.x= 'id', by.y= 'id', all.x=T)
    
  return(full_significant_genes_annotated)
  
}

###STATSLOOP non cd-hit
superstats <- function(data_stats,annotations, column, new_mapping){
    
    full_matrix<- remove_constant_columns(data_stats)
    full_fisher <- fisher_test(full_matrix, new_mapping)
    full_logistic <- logistic_regression(full_matrix,new_mapping)
    full_relationships <- extract_relationships(full_matrix, new_mapping)
    
    
    full_matrix_stats <- merge(full_fisher, full_logistic, by = 'id', all.y = T)
    full_matrix_stats <- merge(full_matrix_stats, full_relationships, by= 'id', all = T)
    full_significant_genes <- extract_significant_genes(full_matrix_stats, 'Full_matrix')
    full_significant_genes_annotated <- merge(full_significant_genes, annotations,by.x= 'id', by.y= 'id', all.x=T)
    
  
  return(full_significant_genes_annotated)
}


###CREATE ABSENCE/PRESENCE PLOT
abs_pres_plot <- function(significant_genes, db){
  if (db == 'MCL'| db =='GCF'){
    col1 <- colnames(significant_genes)[4]
    col2 <- colnames(significant_genes)[5]
  }else{
    col1 <- colnames(significant_genes)[5]
    col2 <- colnames(significant_genes)[6]
  }
  plot <- ggplot(significant_genes[significant_genes$pvalue_adj_fisher < 0.01,], aes_string(x=col1, y= col2,      color='log(pvalue_adj_fisher)' , label = 'id')) +
    geom_point(shape=16) +
    scale_color_viridis(name="log pvalue_adj") +#, trans="reverse") + 
    theme(text = element_text(family="Times"),panel.background = element_rect(fill = "white", colour = "grey50")) + theme_bw() + 
    labs(x = paste0("Presence in " ,col1," bacteria"), y = paste0("Presence in " ,col2," bacteria")) + ggtitle(paste0(db ,' significant genes Absence/Presence fisher test'))
  return(plot)
}

###CREATE VULCANO PLOT
vulcano_plot <- function(data, db){
  to_plot <- data.table(data)
  to_plot$significative_difference <- 'NO'
  
  if (db == 'MCL' | db =='GCF'){
    to_plot[to_plot$pvalue_adj_fisher < 0.01]$significative_difference <- 'YES'
    plot <- ggplot(to_plot, aes(x=log2fold, y=-log10(pvalue_adj_fisher), col=significative_difference, label = id)) + 
      geom_point() + ggtitle(paste0(db ,' significant genes')) +
      theme_minimal() + theme_bw()
  }else{
    to_plot[to_plot$pvalue_adj_log < 0.01]$significative_difference <- 'YES'
    plot <- ggplot(to_plot, aes(x=log2fold, y=-log10(pvalue_adj_log), col=significative_difference, label = id)) + 
      geom_point() + ggtitle(paste0(db ,' significant genes quatity logistic regression')) +
      theme_minimal() + theme_bw()
  }
  return(plot)
}


extract_options <- function(groups, column){
  options <- unique(groups[,column])
  return(options)
}


###Main stats function calling the other functions
stats <- function(input_list, clean_matrix){
  groups <- input_list[[1]]
  column <- input_list[[2]] #Column of mapping file to make the stats with
  db <- input_list[[3]] #data to make the stats
  group1 <- input_list[[4]]
  group2 <- input_list[[5]]
  
  
  
  n_samples <- grep("completeness", colnames(clean_matrix)) - 1
  annotations <- clean_matrix[,c(1,(n_samples+1):ncol(clean_matrix))]
  colnames(annotations)[1] <- 'id'
  
  
  if (db == 'MCL'){
    matrix_stats <- clean_matrix[,1:n_samples] 
    matrix_stats<- remove_constant_columns(matrix_stats)
    if (column == 'Lifestyle'){
      remove_unknowns <- mapping_file[mapping_file$Lifestyle %in% 'Unknown',]
      matrix_stats <- matrix_stats[,-which(names(matrix_stats) %in% remove_unknowns$Sample)]
      
    }
    if (column == 'clusters_hdbscan'){
      remove_unknowns <- groups[groups$clusters_hdbscan %in% 'cluster_Outliers',]
      matrix_stats <- matrix_stats[,-which(names(matrix_stats) %in% remove_unknowns$Sample)]
      
    }
    matrix_annotation <- annotations
    
  }
  if (db == 'KEGG'){
    matrix_stats <- clean_matrix[,c(2:n_samples, grep("keggid", colnames(clean_matrix)))]
    matrix_stats = aggregate(matrix_stats[,1:(n_samples - 1)], by = list(id = matrix_stats[,n_samples]), FUN=sum)
    matrix_stats<- remove_constant_columns(matrix_stats)
    matrix_annotation <- distinct(annotations[,c(grep("keggid", colnames(annotations)), grep("keggid", colnames(annotations)) + 1)])
    
  }
  if (db == 'COG'){
    matrix_stats <- clean_matrix[,c(2:n_samples, grep("cogid", colnames(clean_matrix)))]
    matrix_stats <- aggregate(matrix_stats[,1:(n_samples - 1)], by = list(id = matrix_stats[,n_samples]), FUN=sum)
    matrix_stats<- remove_constant_columns(matrix_stats)
    matrix_annotation <- distinct(annotations[,c(grep("cogid", colnames(annotations)), grep("cogid", colnames(annotations)) + 1)])
  }
  if (db == 'PFAM'){
    matrix_stats <- clean_matrix[,c(2:n_samples, grep("pfamid", colnames(clean_matrix)))]
    matrix_stats <- aggregate(matrix_stats[,1:(n_samples - 1)], by = list(id = matrix_stats[,n_samples]), FUN=sum)
    matrix_stats<- remove_constant_columns(matrix_stats)
    matrix_annotation <- distinct(annotations[,c(grep("pfamid", colnames(annotations)), grep("pfamid", colnames(annotations)) + 1)])
  }
  if (db == 'PROKKA'){
    matrix_stats <- clean_matrix[,c(2:n_samples, grep("prokkadescription", colnames(clean_matrix)))]
    matrix_stats <- matrix_stats[matrix_stats$prokkadescription != 'hypothetical protein',] 
    matrix_stats <- aggregate(matrix_stats[,1:(n_samples - 1)], by = list(id = matrix_stats[,n_samples]), FUN=sum)
    matrix_stats<- remove_constant_columns(matrix_stats)
    matrix_annotation <- distinct(annotations[,c(grep("prokkadescription", colnames(annotations)), grep("prokkadescription", colnames(annotations)))])
    colnames(matrix_annotation)[1] <- 'id'
    colnames(matrix_annotation)[2] <- 'prokkadescription'
    
  }
  if (db == 'DBCAN'){
    matrix_stats <- clean_matrix[,c(2:n_samples, grep("dbcanid", colnames(clean_matrix)))]
    matrix_stats = aggregate(matrix_stats[,1:(n_samples - 1)], by = list(id = matrix_stats[,n_samples]), FUN=sum)
    matrix_stats<- remove_constant_columns(matrix_stats)
    matrix_annotation <- distinct(annotations[,c(grep("dbcanid", colnames(annotations)), grep("dbcanid", colnames(annotations)) + 1)])
  }
  colnames(matrix_annotation)[1] <- 'id'
  
  
  
  new_mapping_file <- data.table(groups[,c('Sample', column)])
  colnames(new_mapping_file)[2]<- 'Lifestyle'
  if (group2 == 'All'){
    new_mapping_file[new_mapping_file$Lifestyle != group1]$Lifestyle <- 'others'
  }else{
    new_mapping_file <-  new_mapping_file[new_mapping_file$Lifestyle == group1|new_mapping_file$Lifestyle == group2]
  }
  
  
  if (db == 'MCL'){
    results <- superstats_cdhit(matrix_stats,matrix_annotation, column, new_mapping_file)
  }else{
    results <- superstats(matrix_stats,matrix_annotation, column, new_mapping_file)
  }
  
  return(results)
  
  
  
  
}

##Function  that fix GCF names and merge with annotations
process_bgc_table <- function(bgc_matrix, bgc_annotation ){
  bgc_matrix$id <- rownames(bgc_matrix)
  bgc_matrix <- bgc_matrix[,c(ncol(bgc_matrix), 1:(ncol(bgc_matrix) - 1))]
  bgc_annotation$GCF.No <- paste0('GCF',bgc_annotation$GCF.No)
  bgc_annotation <- aggregate(  BGC.Class~ GCF.No , data = bgc_annotation, paste, collapse = ",")
  
  M <- merge(bgc_matrix, bgc_annotation, by.x = 'id', by.y = 'GCF.No', all.x = T)
  return(M)
}

###Stats functions for GCF level
stats_gcf <- function(input_list, gcf){
  groups <- input_list[[1]]
  column <- input_list[[2]] #Column of mapping file to make the stats with
  db <- input_list[[3]] #data to make the stats
  group1 <- input_list[[4]]
  group2 <- input_list[[5]]
  
  
  
  matrix_stats <- gcf[,1:grep("BGC.Class", colnames(gcf))] 
  matrix_stats<- remove_constant_columns(matrix_stats)
  matrix_annot <- gcf[,c(1,ncol(gcf))]
  
  
  new_mapping_file <- data.table(groups[,c('Sample', column)])
  colnames(new_mapping_file)[2]<- 'Lifestyle'
  if (group2 == 'All'){
    new_mapping_file[new_mapping_file$Lifestyle != group1]$Lifestyle <- 'others'
  }else{
    new_mapping_file <-  new_mapping_file[new_mapping_file$Lifestyle == group1|new_mapping_file$Lifestyle == group2]
  }
  
  results <- superstats_cdhit(matrix_stats,matrix_annot, column, new_mapping_file)
  return(results)
}

###Function to create interactive heatmap of the BGCs
heatmap_gcf_function <- function(gcf_matrix, mapping_file, column){
  gcf_matrix$BGC.Class <- sapply(strsplit(gcf_matrix$BGC.Class,","), `[`, 1)
  gcf_agg <- aggregate(gcf_matrix[,2:(ncol(gcf_matrix) - 1)], by = list(BGC.Class = gcf_matrix$BGC.Class), FUN= sum)
  
  rownames(gcf_agg) <-  gcf_agg$BGC.Class
  gcf_agg$BGC.Class <- NULL
  gcf_agg <- as.matrix(gcf_agg)
  cols <- glasbey(length(unique(mapping_file[,column])))
  names(cols) <- unique(mapping_file[,column])
  plot <- heatmaply(gcf_agg, col_side_colors  = mapping_file[,column], col_side_palette = cols,showticklabels = c(F,T), height = 800, width = 800)
  return(plot)
}


heatmap_gcf_function2 <- function(gcf_matrix, map, column){
  gcf <- gcf_matrix
  
  rownames(gcf) <- gcf$id
  gcf$id <- NULL
  gcf$BGC.Class <- NULL
  
  gcf <- t(gcf)
  gcf <- gcf[,colSums(gcf)> (nrow(gcf) * 0.01)]
  
  df1 <- data.frame(sample = rownames(gcf), id = 1)
  df2 <- data.frame(sample = mapping_file$Sample, id = 2)
  
  M<- merge(df1, df2, by= 'sample')
  map <- map[map$Sample %in% M$sample,]
  
  gcf <- gcf[map$Sample,]
  
  
  cols <- glasbey(length(unique(map[,column])))
  names(cols) <- unique(map[,column])
  #names(cols)<- c('cluster_1','cluster_2','cluster_3','cluster_4')
  
  plot <- heatmaply(as.matrix(gcf), row_side_colors = map[,column],row_side_palette = cols ,showticklabels = c(F,F), height = 800, width = 1000, ylab = 'Genomes', xlab = 'GCFs')
  return(plot)
}


###When cd-hit stats chosen, this creates a plot showing the COG hierarchy of the significant cd-hit clusters
COG_significative_genes <- function(result, cog, min, db){
  if (db == 'COG'){
    column1 <- 6
    column2 <- 5
    colnames(result)[1] <- 'cogid'
  }
  else {
    column1 <- 5
    column2 <- 4
  }
  
  result$mean_difference <- (result[,colnames(result)[column1]] - result[,colnames(result)[column2]])
  result <- result[result$mean_difference > min | result$mean_difference < -min,]
  result <- result[result$pvalue_adj_fisher < 0.01,]
  
  result <- data.table(result)
  result[is.na(result$cogid)]$cogid <- 'Unknown'
  result <- as.data.frame(result)
  
  M <- merge(result, cog_extended_annotations, by.x = 'cogid', by.y = 'id', all.x = T)
  M <- data.table(M)
  M[is.na(M$processes)]$processes <- 'S'
  M$new <- paste0('UP in ', colnames(M)[6])
  M[M$log2fold >0]$new <- paste0('DOWN in ', colnames(M)[6])
  M <- as.data.frame(M)
  
  df <- as.data.frame(table(M$new,M$processes))
  
  
  M2 <- M
  M2$pvalue_adj_fisher<- -log10(M2$pvalue_adj_fisher)
  
  op <- unique(M$new)
  
  colors <- as.character(alphabet2(n=2))
  color1 <-colors[1]
  color2 <- colors[2]
  fig <- M2 %>%
    plot_ly(type = 'violin')
  fig <- fig %>%
    add_trace(
      x = ~processes[M2$new == op[1]],
      y = ~mean_difference[M2$new == op[1]],
      
      legendgroup = paste0('DOWN in ', colnames(M2)[6]),
      scalegroup = paste0('DOWN in ', colnames(M2)[6]),
      name = paste0('DOWN in ', colnames(M2)[6]),
      
      box = list(
        visible = T
      ),
      meanline = list(
        visible = T
      ),
      color = I('grey')
    )
  fig <- fig %>%
    add_trace(
      x = ~processes[M2$new == op[2]],
      y = ~mean_difference[M2$new == op[2]],
      legendgroup= paste0('UP in ', colnames(M2)[6]),
      scalegroup = paste0('UP in ', colnames(M2)[6]),
      name = paste0('UP in ', colnames(M2)[6]),
      
      box = list(
        visible = T
      ),
      meanline = list(
        visible = T
      ),
      color = I('black')
    )
  
  
  
  fig <- fig %>%
    layout(
      yaxis = list(categoryorder = "array",categoryarray = ~processes,
        zeroline = F, title = "difference in group means"
      )
      
    )
  
  
  
  M2$processes <- factor(M2$processes, levels = unique(M2$processes))
  
  df$Var2 <- factor(df$Var2, levels = unique(M2$processes))
  
  fig2 <- df %>% plot_ly(
    x = ~Var2[df$Var1 == op[1]],
    y = ~Freq[df$Var1 == op[1]],
    
    name = paste0('DOWN in ', colnames(M2)[6]),
    type = "bar",
    color = I('grey')
    
  )
  
  fig2 <- fig2 %>% add_trace(
    x = ~Var2[df$Var1 == op[2]],
    y = ~Freq[df$Var1 == op[2]],
    
    name = paste0('UP in ', colnames(M2)[6]),
    type = "bar",
    color = I('black')
  )
  ax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
  )
  fig2 <- fig2 %>% layout(yaxis = list(title = 'Count'), barmode = 'group', showlegend = FALSE)
  
  
  
  final <- subplot(style(fig2, showlegend = F),
                   style(fig, showlegend = T),
                   nrows = 2, heights = c(0.34,0.66),titleY = TRUE, titleX = F, shareX = T)
  return(final)
}


###Creates the abs/pres plot with all the cd-hit gene clusters
super_abs_pres_plot <- function(mtrx, groups, group1,group2, column){
  new_mapping_file <- data.table(groups[,c('Sample', column)])
  colnames(new_mapping_file)[2]<- 'Lifestyle'
  if (group2 == 'All'){
    new_mapping_file[new_mapping_file$Lifestyle != group1]$Lifestyle <- 'others'
    group2 <- 'others'
  }else{
    new_mapping_file <-  new_mapping_file[new_mapping_file$Lifestyle == group1|new_mapping_file$Lifestyle == group2]
  }
  
  
  
  
  samples_1 <- new_mapping_file[new_mapping_file$Lifestyle == group1]$Sample
  samples_2 <- new_mapping_file[new_mapping_file$Lifestyle == group2]$Sample
  
  #mtrx[2:ncol(mtrx)] <- apply(mtrx[2:ncol(mtrx)], 2, function(x) ifelse(x > 1, 1, x))
  mtrx[2:ncol(mtrx)] <- mtrx[2:ncol(mtrx)] %>% mutate_if(is.numeric, ~1 * (. > 0)) #Faster and less memory
  subset1 <- mtrx[,samples_1]
  subset2 <- mtrx[,samples_2]
  
  group1_means <- rowMeans(subset1[,2:ncol(subset1)])
  group2_means <- rowMeans(subset2[,2:ncol(subset2)])
  
  df <- data.frame(id = mtrx[,1], group1= group1_means, group2 = group2_means)
  colnames(df)[2] <- group1
  colnames(df)[3] <- group2
  
  plot <- ggplot(df, aes_string(x= group1, y = group2)) + geom_hex(bins = 50) + scale_fill_gradient(low="lightblue1",high="darkblue",trans="log10")  + 
    theme_bw() + xlab( paste0('Relative presence in ' , group1)) + ylab( paste0('Relative presence in ', group2))
  plot <- ggplotly(plot)
  return(plot)
}
