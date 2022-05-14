library(data.table)
library(dplyr)
library(ggplot2)
options(bitmapType='cairo')
#shiny::runApp('/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeShiny')

#------------------#
# Define variables #
#------------------#

# Input files
input_dir = '/home/j.aguirreplans/Projects/Scipher/SampleSize/data/data_shiny_app'
plots_dir = '/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/plots_shiny'
topology_results_file = paste(input_dir, 'analysis_topology.csv', sep='/')
ppi_results_file = paste(input_dir, 'analysis_ppi.csv', sep='/')
disease_genes_results_file = paste(input_dir, 'analysis_disease_genes.csv', sep='/')
essential_genes_results_file = paste(input_dir, 'analysis_essential_genes.csv', sep='/')

# Dictionary that maps the name of the variable to a nice name for the plot
parameter2label <- list("nodes"="Nodes", "edges"="Edges", "av_degree"="Av. degree", "av_path_length", "Av. path length", "av_clustering_coef" = "Av. clust. coef.", "num_components" = "Num. of components", "size_lcc" = "Size of the LCC", 
                        "lost_nodes" = "Lost nodes", "lost_edges" = "Lost edges", "gained_nodes" = "Gained nodes", "gained_edges" = "Gained edges",
                        "TP" = "TP", "FP" = "FP", "TN" = "TN", "FN" = "FN", "TPR"="TPR", "FPR"="FPR", "accuracy"="Accuracy", "F1"="F1", "MCC"="MCC", "corr"="Correlation",
                        "disease_genes_in_network" = "Num. disease genes in network", "percent_disease_genes_in_network" = "% disease genes in network", "lcc_size" = "Num. disease genes forming a LCC", "percent_disease_genes_in_lcc" = "% disease genes forming a LCC",
                        "number.of.edges" = "Number of edges", "difference.mean" = "Mean difference between co-expression scores", "score.subset.mean" = "Co-expression score (networks from subsets)", "score.all.samples.mean" = "Co-expression score (networks from all samples)", "distance.mean" = "Mean distance between pairs of proteins (in the PPI network)", "local_clustering_coefficient.mean" = "Mean local clustering coefficient between pairs of proteins (in the PPI network)", "degree_centrality.mean" = "Mean degree between pairs of proteins (in the PPI network)", "betweenness_centrality.mean" = "Mean betweenness centrality between pairs of proteins (in the PPI network)",
                        "num.edges" = "Number of edges", "distances" = "Mean distance between pairs of proteins (in the PPI network)", "difference" = "Mean difference between co-expression scores",
                        "num_nodes" = "Number of nodes", "num_edges" = "Number of significant edges", "num_lcc_nodes" = "Number of nodes in the LCC", "num_lcc_edges" = "Number of edges in the LCC", "lcc_z" = "LCC significance z-score", "lcc_pvalue" = "LCC significance p-value", "log_lcc_pvalue" = "LCC significance log(p-value)",
                        "max_k" = "Maximum k-core", "num_main_core_nodes" = "Number of nodes in the main core", "num_main_core_edges" = "Number of edges in the main core",
                        "num_essential_genes" = "Number of essential genes", "fraction_essential_genes" = "Fraction of essential genes", "num_essential_lcc_nodes" = "Number of essential genes forming a LCC", "fraction_essential_lcc_nodes" = "Essential genes rLCC", 
                        "num_disease_genes" = "Number of disease genes", "fraction_disease_genes" = "Fraction of disease genes", "num_disease_components" = "Number of disease gene components", "num_disease_lcc_nodes" = "Number of disease genes in the LCC", "fraction_disease_lcc_nodes" = "Disease rLCC", "num_disease_lcc_edges" = "Number of disease gene edges in the LCC", "disease_lcc_z" = "Significance z-score of the disease LCC", "disease_lcc_pvalue" = "Significance p-value of the disease LCC", "log_disease_lcc_pvalue" = "Significance log(p-value) of the disease LCC",
                        "overlapindex" = "Overlap index", "jaccardIndex" = "Jaccard index",
                        "num_ppi_nodes" = "Number of PPI nodes", "num_ppi_edges" = "Number of PPI edges", "fraction_ppi_nodes" = "Fraction of PPI nodes", "fraction_ppi_edges" = "Fraction of PPI edges", "num_ppi_main_core_nodes" = "Number of PPI main core nodes", "num_ppi_main_core_edges" = "Number of PPI main core edges", "fraction_ppi_main_core_nodes" = "Fraction of PPI main core nodes", "fraction_ppi_main_core_edges" = "Fraction of PPI main core edges"
)

#-----------#
# Read data #
#-----------#

topology_results_df = fread(topology_results_file)
ppi_results_df = fread(ppi_results_file)
disease_genes_results_df = fread(disease_genes_results_file)
essential_genes_results_df = fread(essential_genes_results_file) %>% rename("num_essential_components" = "num_components", "num_essential_lcc_nodes" = "num_lcc_nodes", "num_essential_lcc_edges" = "num_lcc_edges", "essential_lcc_z" = "lcc_z", "essential_lcc_pvalue" = "lcc_pvalue")
results_df = inner_join(topology_results_df, ppi_results_df, by = c("method", "dataset", "type_dataset", "size", "rep", "type_correlation", "threshold")) %>% inner_join(disease_genes_results_df, by = c("method", "dataset", "type_dataset", "size", "rep", "type_correlation", "threshold")) %>% inner_join(essential_genes_results_df, by = c("method", "dataset", "type_dataset", "size", "rep", "type_correlation", "threshold"))
results_df$type_dataset = tolower(results_df$type_dataset)

#------------------------#
# Send information to UI #
#------------------------#

server <- function(input, output, session) {

  #----------------------------------------------------------#
  # Change information of the UI depending on the selections #
  #----------------------------------------------------------#
  
  observeEvent(input$type_analysis, {
    
    freezeReactiveValue(input, "boxplot_parameter") # Freeze boxplot_parameter until updateSelectInput is completed
    if (input$type_analysis == "topology"){
      updateSelectInput(session = session, 
                        inputId = "boxplot_parameter", 
                        label = "Parameter:", 
                        choices = c("Num. edges" = "num_edges",
                                    "Num. nodes" = "num_nodes",
                                    "Av. degree" = "av_degree",
                                    "Av. path length" = "av_path_length",
                                    "Av. clust. coef." = "av_clustering_coef",
                                    "Num. components" = "num_components",
                                    "Nodes of the LCC" = "num_lcc_nodes",
                                    "Edges of the LCC" = "num_lcc_edges",
                                    "LCC z-score" = "lcc_z",
                                    "LCC p-value" = "lcc_pvalue",
                                    "LCC log(p-value)" = "log_lcc_pvalue",
                                    "Maximum k-core" = "max_k",
                                    "Num. nodes in main core" = "num_main_core_nodes",
                                    "Num. edges in main core" = "num_main_core_edges",
                                    "Repetitions" = "rep"),
                        selected = "num_edges")
    } else if (input$type_analysis == "ppi"){
      updateSelectInput(session = session, 
                        inputId = "boxplot_parameter", 
                        label = "Parameter:", 
                        choices = c('Num. PPI nodes' = 'num_ppi_nodes',
                                    'Num. PPI edges' = 'num_ppi_edges',
                                    'Fraction PPI nodes' = 'fraction_ppi_nodes',
                                    'Fraction PPI edges' = 'fraction_ppi_edges',
                                    'Num. PPI main core nodes' = 'num_ppi_main_core_nodes',
                                    'Num. PPI main core edges' = 'num_ppi_main_core_edges',
                                    'Fraction PPI main core nodes' = 'fraction_ppi_main_core_nodes',
                                    'Fraction PPI main core edges' = 'fraction_ppi_main_core_edges'),
                        selected = "fraction_ppi_edges")
    } else if (input$type_analysis == "disease_genes"){
      updateSelectInput(session = session, 
                        inputId = "boxplot_parameter", 
                        label = "Parameter:", 
                        choices = c('Num. disease genes in network' = 'num_disease_genes',
                                    'Fraction disease genes in network' = 'fraction_disease_genes',
                                    'Num. components of disease genes' = 'num_disease_components',
                                    'Num. disease genes forming a LCC' = 'num_disease_lcc_nodes',
                                    'Disease rLCC' = 'fraction_disease_lcc_nodes',
                                    'Num. edges of the disease gene LCC' = 'num_disease_lcc_edges',
                                    'Disease gene LCC z-score' = 'disease_lcc_z',
                                    'Disease gene LCC p-value' = 'disease_lcc_pvalue',
                                    'Disease gene LCC log(p-value)' = 'log_disease_lcc_pvalue'),
                        selected = "fraction_disease_genes")
    } else if (input$type_analysis == "essential_genes"){
      updateSelectInput(session = session, 
                        inputId = "boxplot_parameter", 
                        label = "Parameter:", 
                        choices = c('Num. essential nodes' = 'num_essential_genes',
                                    'Fraction of essential genes' = 'fraction_essential_genes',
                                    'Num. components' = 'num_essential_components',
                                    'Num. nodes of the LCC' = 'num_essential_lcc_nodes',
                                    'Essential genes rLCC' = 'fraction_essential_lcc_nodes',
                                    'Num. edges of the LCC' = 'num_essential_lcc_edges',
                                    'LCC z-score' = 'essential_lcc_z',
                                    'LCC p-value' = 'essential_lcc_pvalue'),
                        selected = "fraction_essential_genes")
    }
  })

  output$topologyBoxPlot <- renderPlot({
    
    #---------------------------#
    # Parse selected parameters #
    #---------------------------#
    
    # Parse types of dataset (so that if they are not visible, they are not parsed)
    type_datasets_gtex = c()
    if ("gtex" %in% c(input$topology_dataset)){
      topology_gtex_tissues <- c(unlist(strsplit(renderText(input$topology_gtex_tissues, sep='___')(), split='___')))
      topology_gtex_tissues <- gsub(' ', '.', gsub(' - ', '.', gsub('[\\(\\)]', '', topology_gtex_tissues)))
      topology_gtex_sex <- c(unlist(strsplit(renderText(input$topology_gtex_sex, sep='___')(), split='___')))
      type_datasets_gtex = apply(expand.grid(topology_gtex_tissues, topology_gtex_sex), 1, paste, collapse="_")
      type_datasets_gtex <- gsub('_both', '', type_datasets_gtex)
    }
    type_datasets_scipher = c()
    if ("scipher" %in% c(input$topology_dataset)){
      type_datasets_scipher=c(input$topology_type_dataset_scipher)
    }
    type_datasets_tcga = c()
    if ("tcga" %in% c(input$topology_dataset)){
      type_datasets_tcga=c(input$topology_tcga_project)
    }
    type_datasets = c(type_datasets_gtex, type_datasets_scipher, type_datasets_tcga)
      
    #--------------------#
    # Select information #
    #--------------------#

    results_selected_df = results_df %>% filter((dataset %in% c(input$topology_dataset)) & (type_dataset %in% type_datasets) & (method %in% c(input$topology_method)) & (type_correlation %in% c(input$topology_type_correlation)) & (threshold %in% c(input$topology_pvalue_threshold)))
    
    # Select by diseases
    group_vars = c("dataset", "type_dataset", "method", "size", "type_correlation", "threshold")
    if(input$type_analysis == "disease_genes"){
      diseases_selected <- unlist(strsplit(renderText(input$diseases, sep='___')(), split='___'))
      results_selected_df %<>% filter(disease %in% diseases_selected) %>% select(-disease_class) %>% unique()
      group_vars = c(group_vars, "disease")
    } else {
      results_selected_df %<>% select(!(c("num_disease_genes", "num_disease_edges", "fraction_disease_genes", "num_disease_components", "disease", "disease_class", "num_disease_lcc_nodes", "num_disease_lcc_edges", "fraction_disease_lcc_nodes", "disease_lcc_z", "disease_lcc_pvalue"))) %>% unique()
    }
    
    #----------------------------------#
    # Calculate additional information #
    #----------------------------------#
    
    # Calculate mean and standard deviation tendency
    topology_results_selected_by_size_df = results_selected_df %>%
      group_by_at(group_vars) %>%
      summarise_at(vars(input$boxplot_parameter), list(mean=mean, sd=sd)) %>%
      arrange(size)
    
    # Check if the user wants to plot mean or standard deviation
    if(isTRUE(input$topology_sd)){
      combination_metric = "sd"
    } else{
      combination_metric = "mean"
    }
    
    # Check which parameter should be plotted multiple times (if there are multiple checkboxes with multiple elements selected)
    # In this case, the priority is given by the following list (dataset > type_dataset > method...)
    multiple_options_params = c("dataset", "type_dataset", "method", "type_correlation", "threshold", "disease")
    selected_fill_parameter = NULL
    for (multiple_options_param in multiple_options_params){
      if (length(unique(results_selected_df[[multiple_options_param]])) > 1){
        selected_fill_parameter = multiple_options_param
        break
      }
    }

    #----------------------------------#
    # Define main elements of the plot #
    #----------------------------------#
    
    if (is.null(selected_fill_parameter)){
      # Plot without fill
      topology_results_plot = ggplot(results_selected_df, aes(x=size, y=results_selected_df[[input$boxplot_parameter]])) + 
        guides(fill=guide_legend(title="")) +
        geom_point(alpha=0.5, size=3, col=2, fill=2) +
        geom_line(data = topology_results_selected_by_size_df,
                  aes(x = size, y = get(combination_metric), group=1),
                  col=2, lwd=1)
    } else {
      # Plot with fill
      topology_results_plot = ggplot(results_selected_df, aes(x=size, y=results_selected_df[[input$boxplot_parameter]], fill=get(selected_fill_parameter), group=get(selected_fill_parameter), col=get(selected_fill_parameter))) + 
        geom_point(alpha=0.5, size=3) +
        geom_line(data = topology_results_selected_by_size_df,
                  aes(x = size, y = get(combination_metric)),
                  lwd=1)
    }
    
    #----------------#
    # Customize plot #
    #----------------#
    
    # Define labels
    label_x = "Number of samples"
    label_y = parameter2label[[input$boxplot_parameter]]
    
    # Plot standard deviation (if requested by the user)
    if(isTruthy(input$topology_sd)){
      topology_results_plot$layers[[1]] = NULL # Remove geom_point()
      label_y= paste(label_y, " (SD)", sep="") # Change y axis label
    }
    
    # Transform to log scale (if requested by the user)
    if (isTruthy(input$topology_log_x)){
      label_x= paste(label_x, " (Ln)", sep="")
      topology_results_plot = topology_results_plot + 
        scale_x_continuous(trans = scales::log_trans())
    }
    if (isTruthy(input$topology_log_y)){
      label_y= paste(label_y, " (Ln)", sep="")
      topology_results_plot = topology_results_plot + 
        scale_y_continuous(trans = scales::log_trans())
    }
    
    # Customize axes/labels/grid
    topology_results_plot = topology_results_plot + 
      theme_linedraw() +
      xlab(label_x) +
      ylab(label_y) +
      theme(plot.title =  element_text(size = 17, face="bold"), axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 15), legend.text = element_text(size = 14), legend.position="bottom")
    
    topology_results_plot
    
  })

}

