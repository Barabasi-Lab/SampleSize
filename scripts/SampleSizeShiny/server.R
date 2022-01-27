library(data.table)
library(dplyr)
library(ggplot2)
options(bitmapType='cairo')
#renv::restore("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeR")

# Read input data
input_dir <- '/home/j.aguirreplans/Projects/Scipher/SampleSize/data/data_shiny_app'
topology_file = paste(input_dir, 'analysis_topology.csv', sep='/')
disease_genes_file = paste(input_dir, 'analysis_disease_genes.csv', sep='/')
essential_genes_file = paste(input_dir, 'analysis_essential_genes.csv', sep='/')
coexpressed_ppis_file = paste(input_dir, 'coexpressed_ppis.csv', sep='/')
network_similarity_file = paste(input_dir, 'network_similarity_by_size.csv', sep='/')
network_comparison_file = paste(input_dir, 'network_comparison_subsamples_vs_all_samples.csv', sep='/')
difference_file = paste(input_dir, 'difference.csv', sep='/')
scores_file = paste(input_dir, 'scores.csv', sep='/')
threshold_results_file = paste(input_dir, 'threshold.csv', sep='/')
disease_results_file = paste(input_dir, 'disease.csv', sep='/')

topology_df = fread(topology_file)
disease_genes_df = fread(disease_genes_file)
essential_genes_df = fread(essential_genes_file)
coexpressed_ppis_df = fread(coexpressed_ppis_file)
network_similarity_df = fread(network_similarity_file)
network_comparison_df = fread(network_comparison_file)
difference_df = fread(difference_file)
scores_df = fread(scores_file)
threshold_results_df = fread(threshold_results_file)
disease_results_df = fread(disease_results_file)

topology_df$sex[is.na(topology_df$sex)] = "both"
disease_genes_df$sex[is.na(disease_genes_df$sex)] = "both"
essential_genes_df$sex[is.na(essential_genes_df$sex)] = "both"
network_similarity_df$sex[is.na(network_similarity_df$sex)] = "both"
network_comparison_df$sex[is.na(network_comparison_df$sex)] = "both"
topology_df$lcc_pvalue[topology_df$lcc_pvalue == 0] = 0.0000001
topology_df$log_lcc_pvalue = log(topology_df$lcc_pvalue)
disease_genes_df$disease_lcc_pvalue[disease_genes_df$disease_lcc_pvalue == 0] = 0.0000001
disease_genes_df$log_disease_lcc_pvalue = log(disease_genes_df$disease_lcc_pvalue)
essential_genes_df$lcc_pvalue[essential_genes_df$lcc_pvalue == 0] = 0.0000001
essential_genes_df$log_lcc_pvalue = log(essential_genes_df$lcc_pvalue)
network_similarity_df$overlapindex[is.na(network_similarity_df$overlapindex)] = 0.0
network_comparison_df$overlapindex[is.na(network_comparison_df$overlapindex)] = 0.0



# Parameter to label
parameter2label <- list("nodes"="Nodes", "edges"="Edges", "av_degree"="Av. degree", "av_path_length", "Av. path length", "av_clustering_coef" = "Av. clust. coef.", "num_components" = "Num. of components", "size_lcc" = "Size of the LCC", 
                        "lost_nodes" = "Lost nodes", "lost_edges" = "Lost edges", "gained_nodes" = "Gained nodes", "gained_edges" = "Gained edges",
                        "TP" = "TP", "FP" = "FP", "TN" = "TN", "FN" = "FN", "TPR"="TPR", "FPR"="FPR", "accuracy"="Accuracy", "F1"="F1", "MCC"="MCC", "corr"="Correlation",
                        "disease_genes_in_network" = "Num. disease genes in network", "percent_disease_genes_in_network" = "% disease genes in network", "lcc_size" = "Num. disease genes forming a LCC", "percent_disease_genes_in_lcc" = "% disease genes forming a LCC",
                        "number.of.edges" = "Number of edges", "difference.mean" = "Mean difference between co-expression scores", "score.subset.mean" = "Co-expression score (networks from subsets)", "score.all.samples.mean" = "Co-expression score (networks from all samples)", "distance.mean" = "Mean distance between pairs of proteins (in the PPI network)", "local_clustering_coefficient.mean" = "Mean local clustering coefficient between pairs of proteins (in the PPI network)", "degree_centrality.mean" = "Mean degree between pairs of proteins (in the PPI network)", "betweenness_centrality.mean" = "Mean betweenness centrality between pairs of proteins (in the PPI network)",
                        "num.edges" = "Number of edges", "distances" = "Mean distance between pairs of proteins (in the PPI network)", "difference" = "Mean difference between co-expression scores",
                        "num_nodes" = "Number of nodes", "num_edges" = "Number of edges", "num_lcc_nodes" = "Number of nodes in the LCC", "num_lcc_edges" = "Number of edges in the LCC", "lcc_z" = "LCC significance z-score", "lcc_pvalue" = "LCC significance p-value", "log_lcc_pvalue" = "LCC significance log(p-value)",
                        "max_k" = "Maximum k-core", "num_main_core_nodes" = "Number of nodes in the main core", "num_main_core_edges" = "Number of edges in the main core",
                        "num_essential_genes" = "Number of essential genes", "fraction_essential_genes" = "Fraction of essential genes", "num_essential_lcc_nodes" = "Number of essential genes forming a LCC", "fraction_essential_lcc_nodes" = "Fraction of essential genes forming a LCC", 
                        "num_disease_genes" = "Number of disease genes", "fraction_disease_genes" = "Fraction of disease genes", "num_disease_components" = "Number of components formed by genes of the disease", "num_disease_lcc_nodes" = "Number of nodes of the LCC formed by genes of the disease", "fraction_disease_lcc_nodes" = "Fraction of genes in the LCC of the disease", "num_disease_lcc_edges" = "Number of edges of the LCC formed by genes of the disease", "disease_lcc_z" = "Significance z-score of the disease LCC", "disease_lcc_pvalue" = "Significance p-value of the disease LCC", "log_disease_lcc_pvalue" = "Significance log(p-value) of the disease LCC",
                        "overlapindex" = "Overlap index", "jaccardIndex" = "Jaccard index")

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$topologyBoxPlot <- renderPlot({
    
    # To plot using ggplot, have to convert the size vector into character (because if not, it plots a single box wtf)
    topology_df$size = as.character(topology_df$size)
    topology_df$cut = as.character(topology_df$cut)
    topology_df$tissue = tolower(topology_df$tissue)

    get_gtex_tissues <- renderText(input$topology_gtex_tissues, sep='___')
    topology_gtex_tissues <- c(unlist(strsplit(get_gtex_tissues(), split='___')))
    get_gtex_sex <- renderText(input$topology_gtex_sex, sep='___')
    topology_gtex_sex <- c(unlist(strsplit(get_gtex_sex(), split='___')))
    topology_gtex_tissues <- gsub(' ', '.', gsub(' - ', '.', gsub('[\\(\\)]', '', topology_gtex_tissues)))
    
    # Filter by dataset
    if (input$dataset_topology == "gtex"){
      selected_topology_df = topology_df %>% filter((dataset == input$dataset_topology) & (tissue %in% topology_gtex_tissues) & (sex %in% topology_gtex_sex) & (method == input$method_topology))
    } else {
      selected_topology_df = topology_df %>% filter((dataset == input$dataset_topology) & (type_dataset == input$type_scipher_dataset_topology) & (method == input$method_topology))
    }
    
    # Filter samples if there are too many
    sizes = as.numeric(as.vector(selected_topology_df$size[!(selected_topology_df$size == "all")]))
    if (max(sizes) > 500){
      selected_topology_df = selected_topology_df %>% filter(size %in% c(seq(20, max(sizes), 20), "all"))
    }
    selected_topology_df$size <- factor(selected_topology_df$size , levels=c(as.character(sort(as.integer(unique(selected_topology_df$size[!(selected_topology_df$size == "all")])))), "all"))
    
    # Calculate mean
    selected_topology_by_size_df = selected_topology_df %>%
      group_by(size, cut) %>%
      summarise_at(vars(input$parameter_topology), list(mean=mean, median=median, sd=sd, var=var)) %>%
      arrange(size)
    selected_topology_by_size_df$mean.upper = selected_topology_by_size_df$mean + selected_topology_by_size_df$var
    selected_topology_by_size_df$mean.lower = selected_topology_by_size_df$mean - selected_topology_by_size_df$var

    # Plot using ggplot
    if(input$method_topology %in% c("spearman", "pearson", "wto")){
      if(isTruthy(input$pvalue_threshold_topology == 'group_topology_pvalue_threshold')){
        ggplot(selected_topology_df, aes(x=size, y=selected_topology_df[[input$parameter_topology]], fill=cut)) + 
          geom_violin(alpha=0.5) +
          geom_jitter(col=rgb(0.4,0.4,0.8,0.6), pch=16 , cex=1.3, size=1, alpha=0.9) +
          geom_line(data = selected_topology_by_size_df,
                    aes(x = size, y = mean, group=cut, col=cut),
                    lwd=1) +
          theme_linedraw() +
          labs(x="Sample size", y = parameter2label[[input$parameter_topology]])
      } else {
        selected_topology_df = selected_topology_df %>% filter(cut == input$pvalue_threshold_topology)
        selected_topology_by_size_df = selected_topology_by_size_df %>% filter(cut == input$pvalue_threshold_topology)
        ggplot(selected_topology_df, aes(x=size, y=selected_topology_df[[input$parameter_topology]])) + 
          geom_violin(alpha=0.5) +
          scale_x_discrete(limits = levels(selected_topology_df$size)) +
          geom_jitter(col=rgb(0.4,0.4,0.8,0.6), pch=16 , cex=1.3, size=1, alpha=0.9) +
          geom_line(data = selected_topology_by_size_df,
                    aes(x = size, y = mean, group=1),
                    col=2, lwd=1) +
          theme_linedraw() +
          labs(x="Sample size", y = parameter2label[[input$parameter_topology]])
      }
    } else {
      ggplot(selected_topology_df, aes(x=size, y=selected_topology_df[[input$parameter_topology]])) + 
        geom_violin(alpha=0.5) +
        geom_jitter(col=rgb(0.4,0.4,0.8,0.6), pch=16 , cex=1.3, size=1, alpha=0.9) +
        geom_line(data = selected_topology_by_size_df,
                  aes(x = size, y = mean, group=1),
                  col=2, lwd=1) +
        theme_linedraw() +
        labs(x="Sample size", y = parameter2label[[input$parameter_topology]])
    }

    #plot(selected_network_similarity_df$size, selected_network_similarity_df[[input$parameter_similarity]],col=rgb(0.4,0.4,0.8,0.6),pch=16 , cex=1.3, xlab="Sample Size", ylab=parameter2label[[input$parameter_similarity]], main="Similarity between networks with same sample size") 
    #lines(selected_network_similarity_by_size$size, selected_network_similarity_by_size$mean, col=2, lwd=2 )  
    #polygon(c(selected_network_similarity_by_size$size, rev(selected_network_similarity_by_size$size)), c(selected_network_similarity_by_size$mean.upper, rev(selected_network_similarity_by_size$mean.lower)), col = rgb(0.7,0.7,0.7,0.4) , border = NA)

  })

  
  output$essentialityBoxPlot <- renderPlot({
    
    # To plot using ggplot, have to convert the size vector into character (because if not, it plots a single box wtf)
    essential_genes_df$size = as.character(essential_genes_df$size)
    essential_genes_df$cut = as.character(essential_genes_df$cut)
    essential_genes_df$tissue = tolower(essential_genes_df$tissue)
    
    get_gtex_tissues <- renderText(input$essentiality_gtex_tissues, sep='___')
    essentiality_gtex_tissues <- c(unlist(strsplit(get_gtex_tissues(), split='___')))
    get_gtex_sex <- renderText(input$essentiality_gtex_sex, sep='___')
    essentiality_gtex_sex <- c(unlist(strsplit(get_gtex_sex(), split='___')))
    essentiality_gtex_tissues <- gsub(' ', '.', gsub(' - ', '.', gsub('[\\(\\)]', '', essentiality_gtex_tissues)))

    # Filter by dataset
    if (input$dataset_essentiality == "gtex"){
      selected_essential_genes_df = essential_genes_df %>% filter((dataset == input$dataset_essentiality) & (tissue %in% essentiality_gtex_tissues) & (sex %in% essentiality_gtex_sex) & (method == input$method_essentiality))
    } else {
      selected_essential_genes_df = essential_genes_df %>% filter((dataset == input$dataset_essentiality) & (type_dataset == input$type_scipher_dataset_essentiality) & (method == input$method_essentiality))
    }

    # Filter samples if there are too many
    sizes = as.numeric(as.vector(selected_essential_genes_df$size[!(selected_essential_genes_df$size == "all")]))
    if (max(sizes) > 500){
      selected_essential_genes_df = selected_essential_genes_df %>% filter(size %in% c(seq(20, max(sizes), 20), "all"))
    }
    selected_essential_genes_df$size <- factor(selected_essential_genes_df$size , levels=c(as.character(sort(as.integer(unique(selected_essential_genes_df$size[!(selected_essential_genes_df$size == "all")])))), "all"))
    
    # Calculate mean
    selected_essential_genes_by_size_df = selected_essential_genes_df %>%
      group_by(size, cut) %>%
      summarise_at(vars(input$parameter_essentiality), list(mean=mean, median=median, sd=sd, var=var)) %>%
      arrange(size)
    selected_essential_genes_by_size_df$mean.upper = selected_essential_genes_by_size_df$mean + selected_essential_genes_by_size_df$var
    selected_essential_genes_by_size_df$mean.lower = selected_essential_genes_by_size_df$mean - selected_essential_genes_by_size_df$var

    # Plot using ggplot
    if(input$method_essentiality %in% c("spearman", "pearson", "wto")){
      if(isTruthy(input$pvalue_threshold_essentiality == 'group_essentiality_pvalue_threshold')){
        ggplot(selected_essential_genes_df, aes(x=size, y=selected_essential_genes_df[[input$parameter_essentiality]], fill=cut)) + 
          geom_violin(alpha=0.5) +
          geom_jitter(col=rgb(0.4,0.4,0.8,0.6), pch=16 , cex=1.3, size=1, alpha=0.9) +
          geom_line(data = selected_essential_genes_by_size_df,
                    aes(x = size, y = mean, group=cut, col=cut),
                    lwd=1) +
          theme_linedraw() +
          labs(x="Sample size", y = parameter2label[[input$parameter_essentiality]])
      } else {
        selected_essential_genes_df = selected_essential_genes_df %>% filter(cut == input$pvalue_threshold_essentiality)
        selected_essential_genes_by_size_df = selected_essential_genes_by_size_df %>% filter(cut == input$pvalue_threshold_essentiality)
        ggplot(selected_essential_genes_df, aes(x=size, y=selected_essential_genes_df[[input$parameter_essentiality]])) + 
          scale_x_discrete(limits = levels(selected_essential_genes_df$size)) +
          geom_violin(alpha=0.5) +
          geom_jitter(col=rgb(0.4,0.4,0.8,0.6), pch=16 , cex=1.3, size=1, alpha=0.9) +
          geom_line(data = selected_essential_genes_by_size_df,
                    aes(x = size, y = mean, group=1),
                    col=2, lwd=1) +
          theme_linedraw() +
          labs(x="Sample size", y = parameter2label[[input$parameter_essentiality]])
      }
    } else {
      ggplot(selected_essential_genes_df, aes(x=size, y=selected_essential_genes_df[[input$parameter_essentiality]])) + 
        geom_violin(alpha=0.5) +
        geom_jitter(col=rgb(0.4,0.4,0.8,0.6), pch=16 , cex=1.3, size=1, alpha=0.9) +
        geom_line(data = selected_essential_genes_by_size_df,
                  aes(x = size, y = mean, group=1),
                  col=2, lwd=1) +
        theme_linedraw() +
        labs(x="Sample size", y = parameter2label[[input$parameter_essentiality]])
    }
    
    
  })
  
  
  output$diseaseBoxPlot <- renderPlot({
    
    # To plot using ggplot, have to convert the size vector into character (because if not, it plots a single box wtf)
    disease_genes_df$size = as.character(disease_genes_df$size)
    disease_genes_df$cut = as.character(disease_genes_df$cut)
    disease_genes_df$tissue = tolower(disease_genes_df$tissue)
    
    get_gtex_tissues <- renderText(input$disease_genes_gtex_tissues, sep='___')
    disease_genes_gtex_tissues <- c(unlist(strsplit(get_gtex_tissues(), split='___')))
    get_gtex_sex <- renderText(input$disease_genes_gtex_sex, sep='___')
    disease_genes_gtex_sex <- c(unlist(strsplit(get_gtex_sex(), split='___')))
    disease_genes_gtex_tissues <- gsub(' ', '.', gsub(' - ', '.', gsub('[\\(\\)]', '', disease_genes_gtex_tissues)))
    
    get_immune <- renderText(input$immune, sep='___')
    diseases_immune <- unlist(strsplit(get_immune(), split='___'))
    get_endocrine <- renderText(input$endocrine, sep='___')
    diseases_endocrine <- unlist(strsplit(get_endocrine(), split='___'))
    get_cardiovascular <- renderText(input$cardiovascular, sep='___')
    diseases_cardiovascular <- unlist(strsplit(get_cardiovascular(), split='___'))
    get_nervous <- renderText(input$nervous, sep='___')
    diseases_nervous <- unlist(strsplit(get_nervous(), split='___'))
    get_neoplasms <- renderText(input$neoplasms, sep='___')
    diseases_neoplasms <- unlist(strsplit(get_neoplasms(), split='___'))
    get_digestive <- renderText(input$digestive, sep='___')
    diseases_digestive <- unlist(strsplit(get_digestive(), split='___'))
    get_nutritional <- renderText(input$nutritional, sep='___')
    diseases_nutritional <- unlist(strsplit(get_nutritional(), split='___'))
    selected_diseases <- c(diseases_digestive, diseases_nutritional, diseases_nervous, diseases_neoplasms, diseases_immune, diseases_endocrine, diseases_cardiovascular)
    
    # Filter by dataset
    if (input$dataset_disease_genes == "gtex"){
      selected_disease_genes_df = disease_genes_df %>% filter((dataset == input$dataset_disease_genes) & (tissue %in% disease_genes_gtex_tissues) & (sex %in% disease_genes_gtex_sex) & (disease %in% selected_diseases) & (method == input$method_disease_genes))
    } else {
      selected_disease_genes_df = disease_genes_df %>% filter((dataset == input$dataset_disease_genes) & (type_dataset == input$type_scipher_dataset_disease_genes) & (disease %in% selected_diseases) & (method == input$method_disease_genes))
    }
    
    # Filter samples if there are too many
    sizes = as.numeric(as.vector(selected_disease_genes_df$size[!(selected_disease_genes_df$size == "all")]))
    if (max(sizes) > 500){
      selected_disease_genes_df = selected_disease_genes_df %>% filter(size %in% c(seq(20, max(sizes), 20), "all"))
    }
    selected_disease_genes_df$size <- factor(selected_disease_genes_df$size , levels=c(as.character(sort(as.integer(unique(selected_disease_genes_df$size[!(selected_disease_genes_df$size == "all")])))), "all"))
    
    # Calculate mean
    selected_disease_genes_by_size_df = selected_disease_genes_df %>%
      group_by(size, disease, cut) %>%
      summarise_at(vars(input$parameter_disease_genes), list(mean=mean, median=median, sd=sd, var=var)) %>%
      arrange(size)
    selected_disease_genes_by_size_df$mean.upper = selected_disease_genes_by_size_df$mean + selected_disease_genes_by_size_df$var
    selected_disease_genes_by_size_df$mean.lower = selected_disease_genes_by_size_df$mean - selected_disease_genes_by_size_df$var
    
    # Plot using ggplot
    #print(input$group_by)
    if((isTruthy(input$group_by == 'group_disease_pvalue')) & (input$method_disease_genes %in% c("spearman", "pearson", "wto"))){
      selected_disease_genes_df <- selected_disease_genes_df[(selected_disease_genes_df$dataset == input$dataset_disease_genes),]
      ggplot(selected_disease_genes_df, aes(x=size, y=selected_disease_genes_df[[input$parameter_disease_genes]], fill=cut)) + 
        geom_violin(alpha=0.5) +
        geom_jitter(col=rgb(0.4,0.4,0.8,0.6), pch=16 , cex=1.3, size=1, alpha=0.9) +
        geom_line(data = selected_disease_genes_by_size_df,
                  aes(x = size, y = mean, group=cut, col=cut),
                  lwd=1) +
        theme_linedraw() +
        labs(x="Sample size", y = parameter2label[[input$parameter_disease_genes]])
    } else if(isTruthy(input$group_by == 'group_disease')){
      selected_disease_genes_df <- selected_disease_genes_df[(selected_disease_genes_df$dataset == input$dataset_disease_genes) & (selected_disease_genes_df$cut == input$pvalue_threshold_disease_genes),]
      selected_disease_genes_by_size_df = selected_disease_genes_by_size_df %>% filter(cut == input$pvalue_threshold_disease_genes)
      ggplot(selected_disease_genes_df, aes(x=size, y=selected_disease_genes_df[[input$parameter_disease_genes]], fill=disease)) + 
        geom_violin(alpha=0.5) +
        geom_jitter(col=rgb(0.4,0.4,0.8,0.6), pch=16 , cex=1.3, size=1, alpha=0.9) +
        geom_line(data = selected_disease_genes_by_size_df,
                  aes(x = size, y = mean, group=disease, col=disease),
                  lwd=1) +
        theme_linedraw() +
        labs(x="Sample size", y = parameter2label[[input$parameter_disease_genes]])
    } else if(isTruthy(input$group_by == 'group_disease_class')){
      selected_disease_genes_df <- selected_disease_genes_df[(selected_disease_genes_df$dataset == input$dataset_disease_genes) & (selected_disease_genes_df$cut == input$pvalue_threshold_disease_genes),]
      selected_disease_genes_by_size_df = selected_disease_genes_by_size_df %>% filter(cut == input$pvalue_threshold_disease_genes)
      ggplot(selected_disease_genes_df, aes(x=size, y=selected_disease_genes_df[[input$parameter_disease_genes]], fill=disease_class)) + 
        geom_violin(alpha=0.5) +
        geom_jitter(col=rgb(0.4,0.4,0.8,0.6), pch=16 , cex=1.3, size=1, alpha=0.9) +
        geom_line(data = selected_disease_genes_by_size_df,
                  aes(x = size, y = mean, group=disease_class, col=disease_class),
                  lwd=1) +
        theme_linedraw() +
        labs(x="Sample size", y = parameter2label[[input$parameter_disease_genes]])
    } else {
      selected_disease_genes_df <- selected_disease_genes_df[(selected_disease_genes_df$dataset == input$dataset_disease_genes) & (selected_disease_genes_df$cut == input$pvalue_threshold_disease_genes),]
      selected_disease_genes_by_size_df = selected_disease_genes_by_size_df %>% filter(cut == input$pvalue_threshold_disease_genes)
      ggplot(selected_disease_genes_df, aes(x=size, y=selected_disease_genes_df[[input$parameter_disease_genes]])) + 
        geom_violin(alpha=0.5) +
        geom_jitter(col=rgb(0.4,0.4,0.8,0.6), pch=16 , cex=1.3, size=1, alpha=0.9) +
        geom_line(data = selected_disease_genes_by_size_df,
                  aes(x = size, y = mean, group=1),
                  col=2, lwd=1) +
        theme_linedraw() +
        labs(x="Sample size", y = parameter2label[[input$parameter_disease_genes]])
    }
    
  })
  
  
  
  output$diseaseByScoreBoxPlot <- renderPlot({
    
    # To plot using ggplot, have to convert the size vector into character (because if not, it plots a single box wtf)
    disease_results_df$size = factor(disease_results_df$size , levels=as.character(as.numeric(unique(disease_results_df$size))))
    disease_results_df$tissue = tolower(disease_results_df$tissue)
    
    get_gtex_tissues <- renderText(input$disease_genes_by_score_gtex_tissues, sep='___')
    disease_genes_by_score_gtex_tissues <- c(unlist(strsplit(get_gtex_tissues(), split='___')))
    get_gtex_sex <- renderText(input$disease_genes_by_score_gtex_sex, sep='___')
    disease_genes_by_score_gtex_sex <- c(unlist(strsplit(get_gtex_sex(), split='___')))
    disease_genes_by_score_gtex_tissues <- gsub(' ', '.', gsub(' - ', '.', gsub('[\\(\\)]', '', disease_genes_by_score_gtex_tissues)))
    
    get_immune <- renderText(input$disease_genes_by_score_immune, sep='___')
    diseases_immune <- unlist(strsplit(get_immune(), split='___'))
    get_endocrine <- renderText(input$disease_genes_by_score_endocrine, sep='___')
    diseases_endocrine <- unlist(strsplit(get_endocrine(), split='___'))
    get_cardiovascular <- renderText(input$disease_genes_by_score_cardiovascular, sep='___')
    diseases_cardiovascular <- unlist(strsplit(get_cardiovascular(), split='___'))
    get_nervous <- renderText(input$disease_genes_by_score_nervous, sep='___')
    diseases_nervous <- unlist(strsplit(get_nervous(), split='___'))
    get_neoplasms <- renderText(input$disease_genes_by_score_neoplasms, sep='___')
    diseases_neoplasms <- unlist(strsplit(get_neoplasms(), split='___'))
    get_digestive <- renderText(input$disease_genes_by_score_digestive, sep='___')
    diseases_digestive <- unlist(strsplit(get_digestive(), split='___'))
    get_nutritional <- renderText(input$disease_genes_by_score_nutritional, sep='___')
    diseases_nutritional <- unlist(strsplit(get_nutritional(), split='___'))
    selected_diseases <- c(diseases_digestive, diseases_nutritional, diseases_nervous, diseases_neoplasms, diseases_immune, diseases_endocrine, diseases_cardiovascular)

    get_gtex_tissues <- renderText(input$disease_genes_by_score_gtex_tissues, sep='___')
    disease_genes_by_score_gtex_tissues <- c(unlist(strsplit(get_gtex_tissues(), split='___')))
    get_gtex_sex <- renderText(input$disease_genes_by_score_gtex_sex, sep='___')
    disease_genes_by_score_gtex_sex <- c(unlist(strsplit(get_gtex_sex(), split='___')))
    
    # Filter by dataset
    if (input$dataset_disease_genes_by_score == "gtex"){
      selected_disease_results_df <- disease_results_df %>% filter((dataset == input$dataset_disease_genes_by_score) & (tissue %in% disease_genes_by_score_gtex_tissues) & (sex %in% disease_genes_by_score_gtex_sex) & (method == input$method_disease_genes_by_score))
    } else {
      selected_disease_results_df <- disease_results_df %>% filter((dataset == input$dataset_disease_genes_by_score) & (type_dataset == input$type_scipher_dataset_disease_genes_by_score) & (method == input$method_disease_genes_by_score))
    }
    
    # Plot using ggplot
    print(input$disease_genes_by_score_group_by)
    if(isTruthy(input$disease_genes_by_score_group_by == 'group_disease_pval')){
      selected_disease_results_df <- selected_disease_results_df[selected_disease_results_df$disease %in% selected_diseases,]
      ggplot(selected_disease_results_df, aes(x=size, y=.data[[input$metric_disease_genes_by_score]], fill=pval)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$metric_disease_genes_by_score]])
    } else if(isTruthy(input$disease_genes_by_score_group_by == 'group_disease_top')){
      selected_disease_results_df <- selected_disease_results_df[(selected_disease_results_df$disease %in% selected_diseases) & (selected_disease_results_df$pval == input$pval_disease_genes),]
      ggplot(selected_disease_results_df, aes(x=size, y=.data[[input$metric_disease_genes_by_score]], fill=top_percent)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$metric_disease_genes_by_score]])
    } else if(isTruthy(input$disease_genes_by_score_group_by == 'group_disease')){
      selected_disease_results_df <- selected_disease_results_df[selected_disease_results_df$disease %in% selected_diseases,]
      print(dataset)
      print(selected_disease_results_df)
      ggplot(selected_disease_results_df, aes(x=size, y=.data[[input$metric_disease_genes_by_score]], fill=disease)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$metric_disease_genes_by_score]])
    } else if(isTruthy(input$disease_genes_by_score_group_by == 'group_disease_class')){
      selected_disease_results_df <- selected_disease_results_df[(selected_disease_results_df$disease %in% selected_diseases) & (selected_disease_results_df$pval == input$pval_disease_genes) & (selected_disease_results_df$top == input$top_disease_genes),]
      ggplot(selected_disease_results_df, aes(x=size, y=.data[[input$metric_disease_genes_by_score]], fill=disease_class)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$metric_disease_genes_by_score]])
    } else {
      selected_disease_results_df <- selected_disease_results_df[(selected_disease_results_df$disease %in% selected_diseases) & (selected_disease_results_df$pval == input$pval_disease_genes) & (selected_disease_results_df$top == input$top_disease_genes),]
      ggplot(selected_disease_results_df, aes(x=size, y=.data[[input$metric_disease_genes_by_score]])) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$metric_disease_genes_by_score]])
    }
    
  })
  
  
  
  output$similarityBoxPlot <- renderPlot({

    # To plot using ggplot, have to convert the size vector into character (because if not, it plots a single box wtf)
    network_similarity_df$size = as.character(network_similarity_df$size)
    network_similarity_df$cut = as.character(network_similarity_df$cut)
    network_similarity_df$tissue = tolower(network_similarity_df$tissue)
    
    get_gtex_tissues <- renderText(input$similarity_gtex_tissues, sep='___')
    similarity_gtex_tissues <- c(unlist(strsplit(get_gtex_tissues(), split='___')))
    get_gtex_sex <- renderText(input$similarity_gtex_sex, sep='___')
    similarity_gtex_sex <- c(unlist(strsplit(get_gtex_sex(), split='___')))
    similarity_gtex_tissues <- gsub(' ', '.', gsub(' - ', '.', gsub('[\\(\\)]', '', similarity_gtex_tissues)))
    
    # Filter by dataset
    if (input$dataset_similarity == "gtex"){
      selected_network_similarity_df = network_similarity_df %>% filter((dataset == input$dataset_similarity) & (tissue %in% similarity_gtex_tissues) & (sex %in% similarity_gtex_sex) & (method == input$method_similarity) & (type_analysis == input$type_analysis_similarity))
    } else {
      selected_network_similarity_df = network_similarity_df %>% filter((dataset == input$dataset_similarity) & (type_dataset == input$type_scipher_dataset_similarity) & (method == input$method_similarity) & (type_analysis == input$type_analysis_similarity))
    }
    selected_network_similarity_df$size = as.numeric(selected_network_similarity_df$size)

    selected_network_similarity_by_size = selected_network_similarity_df %>%
      group_by(size, cut) %>%
      summarise_at(vars(input$parameter_similarity), list(mean=mean, median=median, sd=sd, var=var)) %>%
      arrange(size)
    selected_network_similarity_by_size$mean.upper = selected_network_similarity_by_size$mean + selected_network_similarity_by_size$var
    selected_network_similarity_by_size$mean.lower = selected_network_similarity_by_size$mean - selected_network_similarity_by_size$var
    
    # Figure based on https://www.r-graph-gallery.com/45-confidence-interval-around-polynomial-curve-fitting.html
    plot(selected_network_similarity_df$size, selected_network_similarity_df[[input$parameter_similarity]],col=rgb(0.4,0.4,0.8,0.6),pch=16 , cex=1.3, xlab="Sample Size", ylab=parameter2label[[input$parameter_similarity]], main="Similarity between networks with same sample size") 
    lines(selected_network_similarity_by_size$size, selected_network_similarity_by_size$mean, col=2, lwd=2 )  
    polygon(c(selected_network_similarity_by_size$size, rev(selected_network_similarity_by_size$size)), c(selected_network_similarity_by_size$mean.upper, rev(selected_network_similarity_by_size$mean.lower)), col = rgb(0.7,0.7,0.7,0.4) , border = NA)
    
  })

  
  
  output$comparisonBoxPlot <- renderPlot({
    network_comparison_df = fread(network_comparison_file)
    
    # To plot using ggplot, have to convert the size vector into character (because if not, it plots a single box wtf)
    network_comparison_df$size = as.character(network_comparison_df$size)
    network_comparison_df$cut = as.character(network_comparison_df$cut)
    network_comparison_df$tissue = tolower(network_comparison_df$tissue)
    
    get_gtex_tissues <- renderText(input$similarity_gtex_tissues, sep='___')
    similarity_gtex_tissues <- c(unlist(strsplit(get_gtex_tissues(), split='___')))
    get_gtex_sex <- renderText(input$similarity_gtex_sex, sep='___')
    similarity_gtex_sex <- c(unlist(strsplit(get_gtex_sex(), split='___')))
    similarity_gtex_tissues <- gsub(' ', '.', gsub(' - ', '.', gsub('[\\(\\)]', '', similarity_gtex_tissues)))
    
    # Filter by dataset
    if (input$dataset_similarity == "gtex"){
      selected_network_comparison_df = network_comparison_df %>% filter((dataset == input$dataset_similarity) & (tissue %in% similarity_gtex_tissues) & (sex %in% similarity_gtex_sex) & (method == input$method_similarity) & (type_analysis == input$type_analysis_similarity))
    } else {
      selected_network_comparison_df = network_comparison_df %>% filter((dataset == input$dataset_similarity) & (type_dataset == input$type_scipher_dataset_similarity) & (method == input$method_similarity) & (type_analysis == input$type_analysis_similarity))
    }
    selected_network_comparison_df$size = as.numeric(selected_network_comparison_df$size)
    
    selected_network_comparison_all_samples_by_size = selected_network_comparison_df %>%
      group_by(size, cut) %>%
      summarise_at(vars(input$parameter_similarity), list(mean=mean, median=median, sd=sd, var=var)) %>%
      arrange(size)
    selected_network_comparison_all_samples_by_size$mean.upper = selected_network_comparison_all_samples_by_size$mean + selected_network_comparison_all_samples_by_size$var
    selected_network_comparison_all_samples_by_size$mean.lower = selected_network_comparison_all_samples_by_size$mean - selected_network_comparison_all_samples_by_size$var
    
    plot(selected_network_comparison_df$size, selected_network_comparison_df[[input$parameter_similarity]],col=rgb(0.4,0.4,0.8,0.6),pch=16 , cex=1.3, xlab="Sample Size", ylab=parameter2label[[input$parameter_similarity]], main="Similarity between subsample networks and the network with all the samples") 
    lines(selected_network_comparison_all_samples_by_size$size, selected_network_comparison_all_samples_by_size$mean, col=2, lwd=2 )  
    polygon(c(selected_network_comparison_all_samples_by_size$size, rev(selected_network_comparison_all_samples_by_size$size)), c(selected_network_comparison_all_samples_by_size$mean.upper, rev(selected_network_comparison_all_samples_by_size$mean.lower)), col = rgb(0.7,0.7,0.7,0.4) , border = NA)
    
  })
  
  
  
  output$coexpressedPPIsBoxPlot <- renderPlot({
    
    # To plot using ggplot, have to convert the size vector into character (because if not, it plots a single box wtf)
    coexpressed_ppis_df$size = as.character(coexpressed_ppis_df$size)
    coexpressed_ppis_df$tissue = tolower(coexpressed_ppis_df$tissue)

    get_gtex_tissues <- renderText(input$coexpressed_ppis_gtex_tissues, sep='___')
    coexpressed_ppis_gtex_tissues <- c(unlist(strsplit(get_gtex_tissues(), split='___')))
    get_gtex_sex <- renderText(input$coexpressed_ppis_gtex_sex, sep='___')
    coexpressed_ppis_gtex_sex <- c(unlist(strsplit(get_gtex_sex(), split='___')))
    coexpressed_ppis_gtex_tissues <- gsub(' ', '.', gsub(' - ', '.', gsub('[\\(\\)]', '', coexpressed_ppis_gtex_tissues)))
    
    # Filter by dataset
    if (input$dataset_coexpressed_ppis == "gtex"){
      selected_coexpressed_ppis_df <- coexpressed_ppis_df %>% filter((type_prediction == input$parameter_coexpressed_ppis) & (dataset == input$dataset_coexpressed_ppis) & (tissue %in% coexpressed_ppis_gtex_tissues) & (sex %in% coexpressed_ppis_gtex_sex) & (method == input$method_coexpressed_ppis))
    } else {
      selected_coexpressed_ppis_df <- coexpressed_ppis_df %>% filter((type_prediction == input$parameter_coexpressed_ppis) & (dataset == input$dataset_coexpressed_ppis) & (type_dataset == input$type_scipher_dataset_coexpressed_ppis) & (method == input$method_coexpressed_ppis))
    }
    selected_coexpressed_ppis_df$size <- factor(selected_coexpressed_ppis_df$size , levels=c(as.character(sort(as.integer(unique(selected_coexpressed_ppis_df$size[(!(selected_coexpressed_ppis_df$size == "all")) & (!(is.na(selected_coexpressed_ppis_df$size)))])))), "all"))

    # Filter samples if there are too many
    sizes = as.numeric(as.vector(selected_coexpressed_ppis_df$size[(!(selected_coexpressed_ppis_df$size == "all")) & (!(is.na(selected_coexpressed_ppis_df$size)))]))
    if (max(sizes) > 500){
      selected_coexpressed_ppis_df = selected_coexpressed_ppis_df %>% filter(size %in% c(seq(20, max(sizes), 20), "all"))
    }
    
    # Calculate mean
    selected_coexpressed_ppis_by_size_df = selected_coexpressed_ppis_df %>%
      group_by(size, cut) %>%
      summarise_at(vars(input$metric_coexpressed_ppis), list(mean=mean, median=median, sd=sd, var=var)) %>%
      arrange(size)
    selected_coexpressed_ppis_by_size_df$mean.upper = selected_coexpressed_ppis_by_size_df$mean + selected_coexpressed_ppis_by_size_df$var
    selected_coexpressed_ppis_by_size_df$mean.lower = selected_coexpressed_ppis_by_size_df$mean - selected_coexpressed_ppis_by_size_df$var

    ggplot(selected_coexpressed_ppis_df, aes(x=size, y=.data[[input$metric_coexpressed_ppis]])) + 
      geom_violin(alpha=0.5) +
      geom_jitter(col=rgb(0.4,0.4,0.8,0.6), pch=16 , cex=1.3, size=1, alpha=0.9) +
      geom_line(data = selected_coexpressed_ppis_by_size_df,
                aes(x = size, y = mean, group=1),
                col=2, lwd=1) +
      theme_linedraw() +
      labs(x="Sample size", y = parameter2label[[input$metric_coexpressed_ppis]])

  })

  output$stabilityBoxPlot <- renderPlot({
    
    # To plot using ggplot, have to convert the size vector into character (because if not, it plots a single box wtf)
    difference_df$range =  paste(difference_df$left.range, difference_df$right.range, sep="-")
    difference_df$size = factor(difference_df$size , levels=as.character(as.numeric(unique(difference_df$size))))
    difference_df$tissue = tolower(difference_df$tissue)
    
    get_gtex_tissues <- renderText(input$stability_gtex_tissues, sep='___')
    stability_gtex_tissues <- c(unlist(strsplit(get_gtex_tissues(), split='___')))
    get_gtex_sex <- renderText(input$stability_gtex_sex, sep='___')
    stability_gtex_sex <- c(unlist(strsplit(get_gtex_sex(), split='___')))
    stability_gtex_tissues <- gsub(' ', '.', gsub(' - ', '.', gsub('[\\(\\)]', '', stability_gtex_tissues)))
    
    # Filter by dataset
    if (input$dataset_stability == "gtex"){
      selected_difference_df <- difference_df %>% filter((dataset == input$dataset_stability) & (tissue %in% stability_gtex_tissues) & (sex %in% stability_gtex_sex) & (method == input$method_stability))
    } else {
      selected_difference_df <- difference_df %>% filter((dataset == input$dataset_stability) & (type_dataset == input$type_scipher_dataset_stability) & (method == input$method_stability))
    }
    
    if(isTruthy(input$stability_group_by == 'group_type_interaction')){
      selected_difference_df <- selected_difference_df %>% filter((type_interaction != "all") & (range == input$difference_range_stability))
      ggplot(selected_difference_df, aes(x=size, y=.data[[input$metric_stability]], fill=type_interaction)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$metric_stability]]) +
        guides(fill=guide_legend(title="Type of interaction"))
    } else if(isTruthy(input$stability_group_by == 'group_difference_range')){
      selected_difference_df <- selected_difference_df %>% filter((range != "0-Inf") & (type_interaction == input$type_interaction_stability))
      #print(selected_difference_df)
      ggplot(selected_difference_df, aes(x=size, y=.data[[input$metric_stability]], fill=range)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$metric_stability]]) + 
        guides(fill=guide_legend(title="Co-expression difference range"))
    } else {
      selected_difference_df <- selected_difference_df %>% filter((type_interaction == input$type_interaction_stability) & (range == input$difference_range_stability))
      #print(selected_difference_df)
      ggplot(selected_difference_df, aes(x=size, y=.data[[input$metric_stability]])) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$metric_stability]])
    }
    
  })

  
  output$scoresBoxPlot <- renderPlot({
    
    # To plot using ggplot, have to convert the size vector into character (because if not, it plots a single box wtf)
    scores_df$range =  paste(scores_df$left.range, scores_df$right.range, sep="-")
    scores_df$size = factor(scores_df$size , levels=as.character(as.numeric(unique(scores_df$size))))
    scores_df$tissue = tolower(scores_df$tissue)
    
    get_gtex_tissues <- renderText(input$scores_gtex_tissues, sep='___')
    scores_gtex_tissues <- c(unlist(strsplit(get_gtex_tissues(), split='___')))
    get_gtex_sex <- renderText(input$scores_gtex_sex, sep='___')
    scores_gtex_sex <- c(unlist(strsplit(get_gtex_sex(), split='___')))
    scores_gtex_tissues <- gsub(' ', '.', gsub(' - ', '.', gsub('[\\(\\)]', '', scores_gtex_tissues)))
    
    # Filter by dataset
    if (input$dataset_scores == "gtex"){
      selected_scores_df <- scores_df %>% filter((dataset == input$dataset_scores) & (tissue %in% scores_gtex_tissues) & (sex %in% scores_gtex_sex) & (method == input$method_scores))
    } else {
      selected_scores_df <- scores_df %>% filter((dataset == input$dataset_scores) & (type_dataset == input$type_scipher_dataset_scores) & (method == input$method_scores))
    }
    
    if (input$type_interaction_scores == "group"){
      if (input$score_range == "group"){
        selected_scores_df <- selected_scores_df %>% filter(type_interaction != "all")
        ggplot(selected_scores_df, aes(x=size, y=.data[[input$metric_scores]], fill=type_interaction)) + 
          geom_boxplot() +
          labs(x="Sample size", y = parameter2label[[input$metric_scores]]) +
          guides(fill=guide_legend(title="Type of interaction"))
      } else {
        selected_scores_df <- selected_scores_df %>% filter((type_interaction != "all") & (range == input$score_range))
        ggplot(selected_scores_df, aes(x=size, y=.data[[input$metric_scores]], fill=type_interaction)) + 
          geom_boxplot() +
          labs(x="Sample size", y = parameter2label[[input$metric_scores]]) +
          guides(fill=guide_legend(title="Type of interaction"))
      }
    } else {
      if ((input$score_range == 'group') | ('group' %in% input$score_range)){
        selected_scores_df <- selected_scores_df %>% filter(type_interaction == input$type_interaction_scores)
        #print(selected_scores_df)
        ggplot(selected_scores_df, aes(x=size, y=.data[[input$metric_scores]], fill=range)) + 
          geom_boxplot() +
          labs(x="Sample size", y = parameter2label[[input$metric_scores]]) +
          guides(fill=guide_legend(title="Co-expression score range"))
      } else {
        if (length(input$score_range) > 1){
          selected_scores_df <- selected_scores_df %>% filter((type_interaction == input$type_interaction_scores) & (range %in% input$score_range))
          #print(selected_scores_df)
          ggplot(selected_scores_df, aes(x=size, y=.data[[input$metric_scores]], fill=range)) + 
            geom_boxplot() +
            labs(x="Sample size", y = parameter2label[[input$metric_scores]]) +
            guides(fill=guide_legend(title="Co-expression score range"))
        } else {
          selected_scores_df <- selected_scores_df %>% filter((type_interaction == input$type_interaction_scores) & (range == input$score_range))
          #print(selected_scores_df)
          ggplot(selected_scores_df, aes(x=size, y=.data[[input$metric_scores]])) + 
            geom_boxplot() +
            labs(x="Sample size", y = parameter2label[[input$metric_scores]])
        }
      }
    }
    
    
  })
  
  
  output$scoresThresholdsBoxPlot <- renderPlot({
    
    # To plot using ggplot, have to convert the size vector into character (because if not, it plots a single box wtf)
    threshold_results_df$threshold = as.character(threshold_results_df$threshold)
    threshold_results_df$size = factor(threshold_results_df$size , levels=as.character(as.numeric(unique(threshold_results_df$size))))
    threshold_results_df$tissue = tolower(threshold_results_df$tissue)
    #print(threshold_results_df)
    
    get_gtex_tissues <- renderText(input$scores_thresholds_gtex_tissues, sep='___')
    scores_thresholds_gtex_tissues <- c(unlist(strsplit(get_gtex_tissues(), split='___')))
    get_gtex_sex <- renderText(input$scores_thresholds_gtex_sex, sep='___')
    scores_thresholds_gtex_sex <- c(unlist(strsplit(get_gtex_sex(), split='___')))
    scores_thresholds_gtex_tissues <- gsub(' ', '.', gsub(' - ', '.', gsub('[\\(\\)]', '', scores_thresholds_gtex_tissues)))
    
    # Filter by dataset
    if (input$dataset_disease_genes_by_score == "gtex"){
      selected_threshold_results_df <- threshold_results_df %>% filter((dataset == input$dataset_scores_thresholds) & (tissue %in% scores_thresholds_gtex_tissues) & (sex %in% scores_thresholds_gtex_sex) & (method == input$method_scores_thresholds))
    } else {
      selected_threshold_results_df <- threshold_results_df %>% filter((dataset == input$dataset_scores_thresholds) & (type_dataset == input$type_scipher_dataset_scores_thresholds) & (method == input$method_scores_thresholds))
    }
    
    if ((input$score_thresholds == 'group') | ('group' %in% input$score_thresholds)){
      selected_threshold_results_df <- selected_threshold_results_df %>% filter(type_interaction == input$type_interaction_scores_thresholds)
      #print(selected_threshold_results_df)
      ggplot(selected_threshold_results_df, aes(x=size, y=.data[[input$metric_scores_thresholds]], fill=threshold)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$metric_scores_thresholds]])
    } else {
      if (length(input$score_thresholds) > 1){
        selected_threshold_results_df <- selected_threshold_results_df %>% filter((type_interaction == input$type_interaction_scores_thresholds) & (threshold %in% input$score_thresholds))
        #print(selected_threshold_results_df)
        ggplot(selected_threshold_results_df, aes(x=size, y=.data[[input$metric_scores_thresholds]], fill=threshold)) + 
          geom_boxplot() +
          labs(x="Sample size", y = parameter2label[[input$metric_scores_thresholds]])
      } else {
        selected_threshold_results_df <- selected_threshold_results_df %>% filter((type_interaction == input$type_interaction_scores_thresholds) & (threshold == input$score_thresholds))
        #print(selected_threshold_results_df)
        ggplot(selected_threshold_results_df, aes(x=size, y=.data[[input$metric_scores_thresholds]])) + 
          geom_boxplot() +
          labs(x="Sample size", y = parameter2label[[input$metric_scores_thresholds]])
      }
    }
    
  })

  
}
