library(data.table)
library(dplyr)
library(ggplot2)
options(bitmapType='cairo')
#renv::restore("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeR")

# Read data
input_dir <- '/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out'
samples_subjects_file = '/home/j.aguirreplans/Databases/GTEx/v8/GTEx_Annotations_SamplesSubjectsMerged.txt'
samples_df = fread(samples_subjects_file)
tissues <- c("Spleen", "Whole.Blood")
sex_to_id = data.frame(row.names=c("male","female") , val=c(1,2))
rep = 10
topology_df <- setNames(data.frame(matrix(ncol = 10, nrow = 0)), c('size', 'rep', 'nodes', 'edges', 'av_degree', 'av_path_length', 'av_clustering_coef', 'pval', 'top', 'dataset'))
components_df <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c('size', 'rep', 'num_components', 'size_lcc', 'pval', 'top', 'dataset'))
composition_df <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), c('size', 'rep', 'lost_nodes', 'lost_edges', 'gained_nodes', 'gained_edges', 'pval', 'top', 'dataset'))
disease_df <- setNames(data.frame(matrix(ncol = 10, nrow = 0)), c('size', 'rep', 'disease_class', 'total_disease_genes', 'total_disease_genes_in_dataset', 'disease_genes_in_network', 'lcc_size', 'pval', 'top', 'dataset'))
coexpressed_ppis_df <- setNames(data.frame(matrix(ncol = 16, nrow = 0)), c("method", "tissue", "sex", "size", "rep", "TP", "FP", "TN", "FN", "TPR", "FPR", "accuracy", "F1", "MCC", "corr", "type_prediction"))
difference_df = data.frame(matrix(ncol=9,nrow=0, dimnames=list(NULL, c("method", "tissue", "sex", "size", "rep", "left.range", "right.range", "number.of.edges", "type_interaction"))))
scores_df = data.frame(matrix(ncol=16,nrow=0, dimnames=list(NULL, c("method", "tissue", "sex", "size", "rep", "left.range", "right.range", "number.of.edges", "difference.mean", "score.subset.mean", "score.all.samples.mean", "distance.mean", "local_clustering_coefficient.mean", "degree_centrality.mean", "betweenness_centrality.mean", "type_interaction"))))
threshold_results_df = data.frame(matrix(ncol=12,nrow=0, dimnames=list(NULL, c("method", "tissue", "sex", "size", "rep", "threshold", "number.of.edges.subset", "number.of.edges.all.samples", "number.of.edges.overlapped", "number.of.edges.lost", "number.of.edges.gained", "type_interaction"))))
disease_results_df = data.frame(matrix(ncol=17,nrow=0, dimnames=list(NULL, c("method", "tissue", "sex", "size", "rep", "disease","disease_class","disease_genes_in_network", "num.edges", "distances", "local_clustering_coefficient.mean", "degree_centrality.mean", "betweenness_centrality.mean", "score.subset.network","score.all.samples.network","difference","type_interaction"))))

for(dataset in c('scipher_nonresponders', 'scipher_responders', 'gtex')){
  if (dataset == 'scipher_nonresponders'){
    results_name = 'analysis_nonresponders_wto_N100'
  } else if (dataset == 'scipher_responders'){
    results_name = 'analysis_responders_wto_N100'
  } else if (dataset == 'gtex'){
    results_name = 'networks_GTEx'
    for(tissue in tissues){
      #tissue.no.sp.char = (samples_df %>% filter(SMTSD == tissue) %>% select(SMTSD.no.sp.char) %>% unique())[[1]]
      for(sex in rownames(sex_to_id)){
        sex_id = sex_to_id[sex,]
        tissue_sex_samples = samples_df[(samples_df$SMTSD.no.sp.char == tissue) & (samples_df$SEX == sex_id),]$SAMPID
        size_list <- seq(10, length(tissue_sex_samples), 10)
        for(size in size_list){
          for(r in 1:rep){
            for(method in c('wto', 'wgcna', 'aracne')){
              coexpressed_ppis_file <- paste(input_dir, '/', results_name, '/', tissue, '_', sex, '/', 'comparison_coexpression_ppi_', method, '_RNAseq_samples_', tissue, '_', sex, '_size_', size, '_rep_', r, '.txt', sep='') # comparison_coexpression_ppi_wgcna_RNAseq_samples_Spleen_female_size_50_rep_10.txt
              if(file.exists(coexpressed_ppis_file)){
                coexpressed_ppis_individual_df <- fread(coexpressed_ppis_file)
                coexpressed_ppis_df <- rbind(coexpressed_ppis_df, cbind(data.frame(method=rep(method, nrow(coexpressed_ppis_individual_df)), tissue=rep(tissue, nrow(coexpressed_ppis_individual_df)), sex=rep(sex, nrow(coexpressed_ppis_individual_df)), size=rep(size, nrow(coexpressed_ppis_individual_df)), rep=rep(r, nrow(coexpressed_ppis_individual_df))), coexpressed_ppis_individual_df))
              }
              difference_file <- paste(input_dir, '/', results_name, '/', tissue, '_', sex, '/', 'analysis_score_difference_', method, '_RNAseq_samples_', tissue, '_', sex, '_size_', size, '_rep_', r, '.txt', sep='') 
              if(file.exists(difference_file)){
                difference_individual_df <- fread(difference_file)
                difference_df <- rbind(difference_df, cbind(data.frame(method=rep(method, nrow(difference_individual_df)), tissue=rep(tissue, nrow(difference_individual_df)), sex=rep(sex, nrow(difference_individual_df)), size=rep(size, nrow(difference_individual_df)), rep=rep(r, nrow(difference_individual_df))), difference_individual_df))
              }
              scores_file <- paste(input_dir, '/', results_name, '/', tissue, '_', sex, '/', 'analysis_score_ranges_', method, '_RNAseq_samples_', tissue, '_', sex, '_size_', size, '_rep_', r, '.txt', sep='') 
              if(file.exists(scores_file)){
                scores_individual_df <- fread(scores_file)
                scores_df <- rbind(scores_df, cbind(data.frame(method=rep(method, nrow(scores_individual_df)), tissue=rep(tissue, nrow(scores_individual_df)), sex=rep(sex, nrow(scores_individual_df)), size=rep(size, nrow(scores_individual_df)), rep=rep(r, nrow(scores_individual_df))), scores_individual_df))
              }
              threshold_results_file <- paste(input_dir, '/', results_name, '/', tissue, '_', sex, '/', 'analysis_score_thresholds_', method, '_RNAseq_samples_', tissue, '_', sex, '_size_', size, '_rep_', r, '.txt', sep='') 
              if(file.exists(threshold_results_file)){
                threshold_results_individual_df <- fread(threshold_results_file)
                threshold_results_df <- rbind(threshold_results_df, cbind(data.frame(method=rep(method, nrow(threshold_results_individual_df)), tissue=rep(tissue, nrow(threshold_results_individual_df)), sex=rep(sex, nrow(threshold_results_individual_df)), size=rep(size, nrow(threshold_results_individual_df)), rep=rep(r, nrow(threshold_results_individual_df))), threshold_results_individual_df))
              }
              disease_results_file <- paste(input_dir, '/', results_name, '/', tissue, '_', sex, '/', 'analysis_score_diseases_', method, '_RNAseq_samples_', tissue, '_', sex, '_size_', size, '_rep_', r, '.txt', sep='') 
              if(file.exists(disease_results_file)){
                disease_results_individual_df <- fread(disease_results_file) %>% select("disease","disease_class","disease_genes_in_network", "num.edges", "distances", "local_clustering_coefficient.mean", "degree_centrality.mean", "betweenness_centrality.mean", "score.subset.network","score.all.samples.network","difference","type_interaction")
                disease_results_df <- rbind(disease_results_df, cbind(data.frame(method=rep(method, nrow(disease_results_individual_df)), tissue=rep(tissue, nrow(disease_results_individual_df)), sex=rep(sex, nrow(disease_results_individual_df)), size=rep(size, nrow(disease_results_individual_df)), rep=rep(r, nrow(disease_results_individual_df))), disease_results_individual_df))
              }
            }
          }
        }
      }
    }
  }
  for(pval in c(0.01, 0.05)){
    for(top_percent in c(0.1, 0.5, 1, 100)){
      cut_off_dir <- paste(input_dir, '/', results_name, '_pval_', pval, '_top_', top_percent, sep='')
      #print(cut_off_dir)
      for(size in c(10, 20, 30, 40, 50)){
        topology_file <- paste(cut_off_dir, '/topology_size_', size, '.txt', sep='')
        components_file <- paste(cut_off_dir, '/components_size_', size, '.txt', sep='')
        composition_file <- paste(cut_off_dir, '/composition_size_', size, '.txt', sep='')
        disease_file <- paste(cut_off_dir, '/disease_components_size_', size, '.txt', sep='')
        if(file.exists(topology_file)){
          small_topology_df <- fread(topology_file)
          small_topology_df$pval <- pval
          small_topology_df$top <- top_percent
          small_topology_df$dataset <- dataset
          topology_df <- rbind(topology_df, small_topology_df)
        }
        if(file.exists(components_file)){
          small_components_df <- fread(components_file)
          small_components_df$pval <- pval
          small_components_df$top <- top_percent
          small_components_df$dataset <- dataset
          components_df <- rbind(components_df, small_components_df)
        }
        if(file.exists(composition_file)){
          small_composition_df <- fread(composition_file)
          small_composition_lost_df <- small_composition_df[small_composition_df$lost_or_gained == 'lost',][,c('size', 'rep', 'nodes', 'edges')]
          colnames(small_composition_lost_df) <- c('size', 'rep', 'lost_nodes', 'lost_edges')
          small_composition_gained_df <- small_composition_df[small_composition_df$lost_or_gained == 'gained',][,c('nodes', 'edges')]
          colnames(small_composition_gained_df) <- c('gained_nodes', 'gained_edges')
          small_composition2_df <- cbind(small_composition_lost_df, small_composition_gained_df)
          small_composition2_df$pval <- pval
          small_composition2_df$top <- top_percent
          small_composition2_df$dataset <- dataset
          composition_df <- rbind(composition_df, small_composition2_df)
        }
        if(file.exists(disease_file)){
          small_disease_df <- fread(disease_file)
          small_disease_df$pval <- pval
          small_disease_df$top <- top_percent
          small_disease_df$dataset <- dataset
          disease_df <- rbind(disease_df, small_disease_df)
        }
      }
    }
  }
}
topology_df <- cbind(topology_df, components_df[,c('num_components', 'size_lcc')], composition_df[,c('lost_nodes', 'lost_edges', 'gained_nodes', 'gained_edges')])

# Calculate the percentage of disease genes in the network and in the LCC
disease_df$percent_disease_genes_in_network = disease_df$disease_genes_in_network/disease_df$total_disease_genes_in_dataset*100
disease_df$percent_disease_genes_in_lcc = disease_df$lcc_size/disease_df$total_disease_genes_in_dataset*100

# Parameter to label
parameter2label <- list("nodes"="Nodes", "edges"="Edges", "av_degree"="Av. degree", "av_path_length", "Av. path length", "av_clustering_coef" = "Av. clust. coef.", "num_components" = "Num. of components", "size_lcc" = "Size of the LCC", 
                        "lost_nodes" = "Lost nodes", "lost_edges" = "Lost edges", "gained_nodes" = "Gained nodes", "gained_edges" = "Gained edges",
                        "TP" = "TP", "FP" = "FP", "TN" = "TN", "FN" = "FN", "TPR"="TPR", "FPR"="FPR", "accuracy"="Accuracy", "F1"="F1", "MCC"="MCC", "corr"="Correlation",
                        "disease_genes_in_network" = "Num. disease genes in network", "percent_disease_genes_in_network" = "% disease genes in network", "lcc_size" = "Num. disease genes forming a LCC", "percent_disease_genes_in_lcc" = "% disease genes forming a LCC",
                        "number.of.edges" = "Number of edges", "difference.mean" = "Mean difference between co-expression scores", "score.subset.mean" = "Co-expression score (networks from subsets)", "score.all.samples.mean" = "Co-expression score (networks from all samples)", "distance.mean" = "Mean distance between pairs of proteins (in the PPI network)", "local_clustering_coefficient.mean" = "Mean local clustering coefficient between pairs of proteins (in the PPI network)", "degree_centrality.mean" = "Mean degree between pairs of proteins (in the PPI network)", "betweenness_centrality.mean" = "Mean betweenness centrality between pairs of proteins (in the PPI network)",
                        "num.edges" = "Number of edges", "distances" = "Mean distance between pairs of proteins (in the PPI network)", "difference" = "Mean difference between co-expression scores")

# Define server logic required to draw a histogram
server <- function(input, output) {

  # Check the selection of group
  observe({
    if (input$pval_topology == "group_topology_pval") {
      #shinyjs::disable(selector = "[type=radio][id=top_percent][value=group_top]")
      #shinyjs::runjs("$('[type=radio][id=top_percent][value=group_top]').parent().parent().addClass('disabled').css('opacity', 0.4)")
      #shinyjs::runjs("$('input[value=group_top]').parent().attr('disabled', true);")
      #shinyjs::disable(selector = "#top_percent button:eq(1)")
      #print('HI1')
    } else if (input$top_topology == "group_topology_top"){
      #shinyjs::disable(selector = "[type=radio][id=pval][value=group_pval]")
      #shinyjs::runjs("$('[type=radio][id=pval][value=group_pval]').parent().parent().addClass('disabled').css('opacity', 0.4)")
      #shinyjs::disable(selector = "#pval")
      #print('HI2')
    } else {
      #shinyjs::enable(selector = "[id=pval]")
      #shinyjs::enable(selector = "[id=top_percent]")
      #shinyjs::runjs("$('[type=radio][id=pval][value=group]').parent().parent().addClass('disabled').css('opacity', 0.4)")
      #shinyjs::disable(selector = "#pval")
      #shinyjs::disable(selector = "[id=top_percent]")
      #print('HI3')
    }
  })

  output$topologyBoxPlot <- renderPlot({
    
    # To plot using ggplot, have to convert the size vector into character (because if not, it plots a single box wtf)
    topology_df$size = as.character(topology_df$size)
    topology_df$pval = as.character(topology_df$pval)
    topology_df$top = as.character(topology_df$top)
    
    # Plot using ggplot
    if(isTruthy(input$pval_topology == 'group_topology_pval')){
      selected_topology_df <- topology_df[(topology_df$dataset == input$dataset_topology) & (topology_df$top == input$top_topology),]
      ggplot(selected_topology_df, aes(x=size, y=selected_topology_df[[input$parameter_topology]], fill=pval)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$parameter_topology]])
    } else if(isTruthy(input$top == 'group_topology_top')){
      selected_topology_df <- topology_df[(topology_df$dataset == input$dataset_topology) & (topology_df$pval == input$pval_topology),]
      ggplot(selected_topology_df, aes(x=size, y=selected_topology_df[[input$parameter_topology]], fill=top_percent)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$parameter_topology]])
    } else {
      selected_topology_df <- topology_df[(topology_df$dataset == input$dataset_topology) & (topology_df$pval == input$pval_topology) & (topology_df$top == input$top_topology),]
      ggplot(selected_topology_df, aes(x=size, y=selected_topology_df[[input$parameter_topology]])) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$parameter_topology]])
    }

    # Using R function of boxplot
    #boxplot(selected_topology_df[[input$parameter_topology]]~size,data=selected_topology_df,
    #        xlab="Sample size", ylab="Number of edges")
  
  })

  output$coexpressedPPIsBoxPlot <- renderPlot({
    
    # To plot using ggplot, have to convert the size vector into character (because if not, it plots a single box wtf)
    coexpressed_ppis_df$size = as.character(coexpressed_ppis_df$size)
    coexpressed_ppis_df$size <- factor(coexpressed_ppis_df$size , levels=as.character(as.numeric(unique(coexpressed_ppis_df$size))))
    coexpressed_ppis_df$tissue = tolower(coexpressed_ppis_df$tissue)

    get_gtex_tissues <- renderText(input$coexpressed_ppis_gtex_tissues, sep='___')
    coexpressed_ppis_gtex_tissues <- c(unlist(strsplit(get_gtex_tissues(), split='___')))
    get_gtex_sex <- renderText(input$coexpressed_ppis_gtex_sex, sep='___')
    coexpressed_ppis_gtex_sex <- c(unlist(strsplit(get_gtex_sex(), split='___')))

    selected_coexpressed_ppis_df <- coexpressed_ppis_df %>% filter((type_prediction == input$parameter_coexpressed_ppis) & (tissue %in% coexpressed_ppis_gtex_tissues) & (sex %in% coexpressed_ppis_gtex_sex) & (method == input$method_coexpressed_ppis))
    #print(selected_coexpressed_ppis_df)
    ggplot(selected_coexpressed_ppis_df, aes(x=size, y=.data[[input$metric_coexpressed_ppis]])) + 
      geom_boxplot() +
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

    if(isTruthy(input$stability_group_by == 'group_type_interaction')){
      selected_difference_df <- difference_df %>% filter((type_interaction != "all") & (range == input$difference_range_stability) & (tissue %in% stability_gtex_tissues) & (sex %in% stability_gtex_sex) & (method == input$method_stability))
      ggplot(selected_difference_df, aes(x=size, y=.data[[input$metric_stability]], fill=type_interaction)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$metric_stability]]) +
        guides(fill=guide_legend(title="Type of interaction"))
    } else if(isTruthy(input$stability_group_by == 'group_difference_range')){
      selected_difference_df <- difference_df %>% filter((range != "0-Inf") & (type_interaction == input$type_interaction_stability) & (tissue %in% stability_gtex_tissues) & (sex %in% stability_gtex_sex) & (method == input$method_stability))
      #print(selected_difference_df)
      ggplot(selected_difference_df, aes(x=size, y=.data[[input$metric_stability]], fill=range)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$metric_stability]]) + 
        guides(fill=guide_legend(title="Co-expression difference range"))
    } else {
      selected_difference_df <- difference_df %>% filter((type_interaction == input$type_interaction_stability) & (range == input$difference_range_stability) & (tissue %in% stability_gtex_tissues) & (sex %in% stability_gtex_sex) & (method == input$method_stability))
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
    
    if (input$type_interaction_scores == "group"){
      if (input$score_range == "group"){
        selected_scores_df <- scores_df %>% filter((type_interaction != "all") & (tissue %in% scores_gtex_tissues) & (sex %in% scores_gtex_sex) & (method == input$method_scores))
        ggplot(selected_scores_df, aes(x=size, y=.data[[input$metric_scores]], fill=type_interaction)) + 
          geom_boxplot() +
          labs(x="Sample size", y = parameter2label[[input$metric_scores]]) +
          guides(fill=guide_legend(title="Type of interaction"))
      } else {
        selected_scores_df <- scores_df %>% filter((type_interaction != "all") & (range == input$score_range) & (tissue %in% scores_gtex_tissues) & (sex %in% scores_gtex_sex) & (method == input$method_scores))
        ggplot(selected_scores_df, aes(x=size, y=.data[[input$metric_scores]], fill=type_interaction)) + 
          geom_boxplot() +
          labs(x="Sample size", y = parameter2label[[input$metric_scores]]) +
          guides(fill=guide_legend(title="Type of interaction"))
      }
    } else {
      if ((input$score_range == 'group') | ('group' %in% input$score_range)){
        selected_scores_df <- scores_df %>% filter((type_interaction == input$type_interaction_scores) & (tissue %in% scores_gtex_tissues) & (sex %in% scores_gtex_sex) & (method == input$method_scores))
        #print(selected_scores_df)
        ggplot(selected_scores_df, aes(x=size, y=.data[[input$metric_scores]], fill=range)) + 
          geom_boxplot() +
          labs(x="Sample size", y = parameter2label[[input$metric_scores]]) +
          guides(fill=guide_legend(title="Co-expression score range"))
      } else {
        if (length(input$score_range) > 1){
          selected_scores_df <- scores_df %>% filter((type_interaction == input$type_interaction_scores) & (range %in% input$score_range)  & (tissue %in% scores_gtex_tissues) & (sex %in% scores_gtex_sex) & (method == input$method_scores))
          #print(selected_scores_df)
          ggplot(selected_scores_df, aes(x=size, y=.data[[input$metric_scores]], fill=range)) + 
            geom_boxplot() +
            labs(x="Sample size", y = parameter2label[[input$metric_scores]]) +
            guides(fill=guide_legend(title="Co-expression score range"))
        } else {
          selected_scores_df <- scores_df %>% filter((type_interaction == input$type_interaction_scores) & (range == input$score_range) & (tissue %in% scores_gtex_tissues) & (sex %in% scores_gtex_sex) & (method == input$method_scores))
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
    
    if ((input$score_thresholds == 'group') | ('group' %in% input$score_thresholds)){
      selected_threshold_results_df <- threshold_results_df %>% filter((type_interaction == input$type_interaction_scores_thresholds) & (tissue %in% scores_thresholds_gtex_tissues) & (sex %in% scores_thresholds_gtex_sex) & (method == input$method_scores_thresholds))
      #print(selected_threshold_results_df)
      ggplot(selected_threshold_results_df, aes(x=size, y=.data[[input$metric_scores_thresholds]], fill=threshold)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$metric_scores_thresholds]])
    } else {
      if (length(input$score_thresholds) > 1){
        selected_threshold_results_df <- threshold_results_df %>% filter((type_interaction == input$type_interaction_scores_thresholds) & (threshold %in% input$score_thresholds)  & (tissue %in% scores_thresholds_gtex_tissues) & (sex %in% scores_thresholds_gtex_sex) & (method == input$method_scores_thresholds))
        #print(selected_threshold_results_df)
        ggplot(selected_threshold_results_df, aes(x=size, y=.data[[input$metric_scores_thresholds]], fill=threshold)) + 
          geom_boxplot() +
          labs(x="Sample size", y = parameter2label[[input$metric_scores_thresholds]])
      } else {
        selected_threshold_results_df <- threshold_results_df %>% filter((type_interaction == input$type_interaction_scores_thresholds) & (threshold == input$score_thresholds) & (tissue %in% scores_thresholds_gtex_tissues) & (sex %in% scores_thresholds_gtex_sex) & (method == input$method_scores_thresholds))
        #print(selected_threshold_results_df)
        ggplot(selected_threshold_results_df, aes(x=size, y=.data[[input$metric_scores_thresholds]])) + 
          geom_boxplot() +
          labs(x="Sample size", y = parameter2label[[input$metric_scores_thresholds]])
      }
    }
    
  })

  
  output$diseaseByScoreBoxPlot <- renderPlot({
    
    # To plot using ggplot, have to convert the size vector into character (because if not, it plots a single box wtf)
    disease_results_df$size = factor(disease_results_df$size , levels=as.character(as.numeric(unique(disease_results_df$size))))
    disease_results_df$tissue = tolower(disease_results_df$tissue)

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
    print(selected_diseases)
    
    get_gtex_tissues <- renderText(input$disease_genes_by_score_gtex_tissues, sep='___')
    disease_genes_by_score_gtex_tissues <- c(unlist(strsplit(get_gtex_tissues(), split='___')))
    get_gtex_sex <- renderText(input$disease_genes_by_score_gtex_sex, sep='___')
    disease_genes_by_score_gtex_sex <- c(unlist(strsplit(get_gtex_sex(), split='___')))
    
    # Plot using ggplot
    print(input$disease_genes_by_score_group_by)
    if(isTruthy(input$disease_genes_by_score_group_by == 'group_disease_pval')){
      selected_disease_results_df <- disease_results_df[(disease_results_df$disease %in% selected_diseases) & (tissue %in% input$disease_genes_by_score_gtex_tissues) & (sex %in% input$disease_genes_by_score_gtex_sex) & (method == input$method_disease_genes_by_score),]
      ggplot(selected_disease_results_df, aes(x=size, y=.data[[input$metric_disease_genes_by_score]], fill=pval)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$metric_disease_genes_by_score]])
    } else if(isTruthy(input$disease_genes_by_score_group_by == 'group_disease_top')){
      selected_disease_results_df <- disease_results_df[(disease_results_df$disease %in% selected_diseases) & (disease_results_df$pval == input$pval_disease_genes) & (tissue %in% input$disease_genes_by_score_gtex_tissues) & (sex %in% input$disease_genes_by_score_gtex_sex) & (method == input$method_disease_genes_by_score),]
      ggplot(selected_disease_results_df, aes(x=size, y=.data[[input$metric_disease_genes_by_score]], fill=top_percent)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$metric_disease_genes_by_score]])
    } else if(isTruthy(input$disease_genes_by_score_group_by == 'group_disease')){
      selected_disease_results_df <- disease_results_df[(disease_results_df$disease %in% selected_diseases) & (disease_results_df$tissue %in% input$disease_genes_by_score_gtex_tissues) & (disease_results_df$sex %in% input$disease_genes_by_score_gtex_sex) & (disease_results_df$method == input$method_disease_genes_by_score),]
      print(tissue)
      print(selected_disease_results_df)
      ggplot(selected_disease_results_df, aes(x=size, y=.data[[input$metric_disease_genes_by_score]], fill=disease)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$metric_disease_genes_by_score]])
    } else if(isTruthy(input$disease_genes_by_score_group_by == 'group_disease_class')){
      selected_disease_results_df <- disease_results_df[(disease_results_df$disease %in% selected_diseases) & (disease_results_df$pval == input$pval_disease_genes) & (disease_results_df$top == input$top_disease_genes) & (tissue %in% input$disease_genes_by_score_gtex_tissues) & (sex %in% input$disease_genes_by_score_gtex_sex) & (method == input$method_disease_genes_by_score),]
      ggplot(selected_disease_results_df, aes(x=size, y=.data[[input$metric_disease_genes_by_score]], fill=disease_class)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$metric_disease_genes_by_score]])
    } else {
      selected_disease_results_df <- disease_results_df[(disease_results_df$disease %in% selected_diseases) & (disease_results_df$pval == input$pval_disease_genes) & (disease_results_df$top == input$top_disease_genes) & (tissue %in% input$disease_genes_by_score_gtex_tissues) & (sex %in% input$disease_genes_by_score_gtex_sex) & (method == input$method_disease_genes_by_score),]
      ggplot(selected_disease_results_df, aes(x=size, y=.data[[input$metric_disease_genes_by_score]])) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$metric_disease_genes_by_score]])
    }
    
  })
  
  
  output$diseaseBoxPlot <- renderPlot({
    
    # To plot using ggplot, have to convert the size vector into character (because if not, it plots a single box wtf)
    disease_df$size = as.character(disease_df$size)
    disease_df$pval = as.character(disease_df$pval)
    disease_df$top = as.character(disease_df$top)
    
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

    # Plot using ggplot
    #print(input$group_by)
    if(isTruthy(input$group_by == 'group_disease_pval')){
      selected_disease_df <- disease_df[(disease_df$dataset == input$dataset_disease_genes) & (disease_df$disease %in% selected_diseases) & (disease_df$top == input$top_disease_genes),]
      ggplot(selected_disease_df, aes(x=size, y=selected_disease_df[[input$parameter_disease_genes]], fill=pval)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$parameter_disease_genes]])
    } else if(isTruthy(input$group_by == 'group_disease_top')){
      selected_disease_df <- disease_df[(disease_df$dataset == input$dataset_disease_genes) & (disease_df$disease %in% selected_diseases) & (disease_df$pval == input$pval_disease_genes),]
      ggplot(selected_disease_df, aes(x=size, y=selected_disease_df[[input$parameter_disease_genes]], fill=top_percent)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$parameter_disease_genes]])
    } else if(isTruthy(input$group_by == 'group_disease')){
      selected_disease_df <- disease_df[(disease_df$dataset == input$dataset_disease_genes) & (disease_df$disease %in% selected_diseases) & (disease_df$pval == input$pval_disease_genes) & (disease_df$top == input$top_disease_genes),]
      ggplot(selected_disease_df, aes(x=size, y=selected_disease_df[[input$parameter_disease_genes]], fill=disease)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$parameter_disease_genes]])
    } else if(isTruthy(input$group_by == 'group_disease_class')){
      selected_disease_df <- disease_df[(disease_df$dataset == input$dataset_disease_genes) & (disease_df$disease %in% selected_diseases) & (disease_df$pval == input$pval_disease_genes) & (disease_df$top == input$top_disease_genes),]
      ggplot(selected_disease_df, aes(x=size, y=selected_disease_df[[input$parameter_disease_genes]], fill=disease_class)) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$parameter_disease_genes]])
    } else {
      selected_disease_df <- disease_df[(disease_df$dataset == input$dataset_disease_genes) & (disease_df$disease %in% selected_diseases) & (disease_df$pval == input$pval_disease_genes) & (disease_df$top == input$top_disease_genes),]
      ggplot(selected_disease_df, aes(x=size, y=selected_disease_df[[input$parameter_disease_genes]])) + 
        geom_boxplot() +
        labs(x="Sample size", y = parameter2label[[input$parameter_disease_genes]])
    }

  })

}
