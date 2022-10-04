library(data.table)
library(dplyr)
library(ggplot2)
library(stats)
#shiny::runApp('/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeShiny')
#shiny::runApp('/Users/j.aguirreplans/Dropbox (CCNR)/Biology/Quim/Scipher/SampleSize/scripts/SampleSizeShiny')
options(bitmapType='cairo') # Only for RStudio webserver if I want to save a plot

#------------------#
# Define variables #
#------------------#

# Input files
#input_dir = '/home/j.aguirreplans/Projects/Scipher/SampleSize/data/data_shiny_app'
#setwd("/Users/j.aguirreplans/Dropbox (CCNR)/Biology/Quim/Scipher/SampleSize/scripts/SampleSizeShiny")
#setwd("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeShiny")
input_dir = 'data'
topology_results_file = paste(input_dir, 'analysis_topology.csv', sep='/')
ppi_results_file = paste(input_dir, 'analysis_ppi.csv', sep='/')
disease_genes_results_file = paste(input_dir, 'analysis_disease_genes.csv', sep='/')
essential_genes_results_file = paste(input_dir, 'analysis_essential_genes.csv', sep='/')
numbers_complete_graph_file = paste(input_dir, 'dataset_numbers_complete_graph.txt', sep='/')

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
                        "num_disease_genes" = "Number of disease genes", "fraction_disease_genes" = "Fraction of disease genes", "num_disease_components" = "Number of disease gene components", "num_disease_lcc_nodes" = "Number of disease genes in the LCC", "fraction_disease_lcc_nodes" = "Disease rLCC", "num_disease_lcc_edges" = "Number of disease gene edges in the LCC", "disease_lcc_z" = "Significance z-score of the disease LCC", "disease_lcc_pvalue" = "Significance p-value of the disease LCC", "log_disease_lcc_pvalue" = "Disease LCC log p-value (abs)",
                        "overlapindex" = "Overlap index", "jaccardIndex" = "Jaccard index",
                        "num_ppi_nodes" = "Number of PPI nodes", "num_ppi_edges" = "Number of PPI edges", "fraction_ppi_nodes" = "Fraction of PPI nodes", "fraction_ppi_edges" = "Fraction of PPI edges", "num_ppi_main_core_nodes" = "Number of PPI main core nodes", "num_ppi_main_core_edges" = "Number of PPI main core edges", "fraction_ppi_main_core_nodes" = "Fraction of PPI main core nodes", "fraction_ppi_main_core_edges" = "Fraction of PPI main core edges",
                        "dataset" = "Dataset", "type_dataset" = "Dataset", "method" = "Method", "type_correlation" = "Type of correlation", "threshold" = "P-value threshold", "disease" = "Disease", "model" = "Model type"
)

#------------------#
# Define functions #
#------------------#

#'  get_logarithmic_tendency
#'  Method to obtain a logarithmic fit from the data.
#'  @param results_dataframe Dataframe containing the data.
#'  @param y_parameter Parameter of the Y axis.
#'  @param x_parameter Parameter of the X axis.
#'  
get_logarithmic_tendency = function(results_dataframe, y_parameter, x_parameter){
  lm_summary = summary(lm(get(y_parameter)~log(get(x_parameter)), data=results_dataframe))
  used_data=data.frame(y=results_dataframe[[y_parameter]], x=log(results_dataframe[[x_parameter]]))
  return(list(lm_summary=lm_summary, used_data=used_data, a=coef(lm_summary)[2], b=coef(lm_summary)[1], L=NaN, adj.r.squared=lm_summary$adj.r.squared))
}

#'  calculate_stretched_exponential_model_without_L
#'  Method to obtain an stretched exponential from the data without L.
#'  @param results_dataframe Dataframe containing the data.
#'  @param y_parameter Parameter of the Y axis.
#'  @param x_parameter Parameter of the X axis.
#'  
calculate_stretched_exponential_model_without_L = function(results_dataframe, y_parameter, x_parameter){
  # Calculate mean of repetitions from same sample size
  results_dataframe_mean = results_dataframe %>% 
    arrange(get(x_parameter), rep) %>%
    group_by(get(x_parameter)) %>%
    summarise_at(vars(all_of(y_parameter)), list(mean=mean, median=median, sd=sd)) %>%
    rename(!!x_parameter := "get(x_parameter)") %>%
    ungroup()
  # Calculate difference between consecutive sizes and y_parameter
  diff_n = diff(results_dataframe_mean[[x_parameter]])
  diff_s = diff(results_dataframe_mean$mean)
  # Calculate gradient = ratio between difference of y_parameter and size
  grad = 1/(diff_n/diff_s)
  # Calculate normalized gradient dividing it by size (starting from 2nd position)
  frac_grad = grad/results_dataframe_mean$mean[2:length(results_dataframe_mean$mean)]
  # Check which factions are negative so that we ignore them. 
  # Why? Because some consecutive differences might give negative values! 
  # This does not make sense in theory, and as the number of situations like this is very small, we ignore it.
  boo=frac_grad>0
  # Calculate polynomial equation = logarithm of the gradient vs. logarithm of sample size (starting from 2nd position) (removing negatives)
  lm_summary = summary(lm(log(frac_grad[boo])~log(results_dataframe_mean[[x_parameter]][2:length(results_dataframe_mean[[x_parameter]])][boo])))
  used_data = data.frame(y=log(frac_grad[boo]), x=log(results_dataframe_mean[[x_parameter]][2:length(results_dataframe_mean[[x_parameter]])][boo]))
  return(list(lm_summary=lm_summary, used_data=used_data, a=coef(lm_summary)[2], b=coef(lm_summary)[1], L=NaN, adj.r.squared=lm_summary$adj.r.squared))
}

#'  calculate_predictions_using_stretched_exponential_model_without_L
#'  Formula to calculate stretched exponential of significant interactions from a list of sample sizes without using L.
#'  This formula does not require the use of L
#'  @param x List of sample sizes.
#'  @param a Slope coefficient.
#'  @param b Intercept coefficient.
#'  
calculate_predictions_using_stretched_exponential_model_without_L = function(x, a, b){
  y = (exp( (exp(b) * x**(1 + a) ) / (1 + a) ))
  y=y/max(y) # Normalize y vector
  return(y)
}

#'  calculate_stretched_exponential_model_by_linear_fit
#'  Formula to find linear fit between log(ln(L)-ln(s)) and log(N).
#'  @param results_dataframe Dataframe containing the data.
#'  @param y_parameter Parameter of the Y axis.
#'  @param x_parameter Parameter of the X axis.
#'  @param L Value of L.
#'  
calculate_stretched_exponential_model_by_linear_fit = function(results_dataframe, y_parameter, x_parameter, L){
  # Calculate mean of repetitions from same sample size
  results_dataframe_mean = results_dataframe %>% 
    arrange(get(x_parameter), rep) %>%
    group_by(get(x_parameter)) %>%
    summarise_at(vars(all_of(y_parameter)), list(mean=mean, median=median, sd=sd)) %>%
    rename(!!x_parameter := "get(x_parameter)") %>%
    ungroup()
  # Define s and N
  s = results_dataframe_mean$mean / max(results_dataframe_mean$mean)
  N = results_dataframe_mean$size
  # Calculate ln(L)-ln(s)
  s_rec = log(L) - log(s)
  # Calculate linear regression between ln(L)-ln(s) and ln(N)
  lm_summary = summary(lm(log(s_rec)~log(N)))
  # Calculate alpha and b
  a = 1-coef(lm_summary)[2]
  b = exp((coef(lm_summary)[1] + log(a-1)))
  # This is the data used to obtain the final model
  used_data = data.frame(y=log(s_rec), x=log(N))
  #return(list(lm_summary=lm_summary, used_data=used_data, a=coef(lm_summary)[2], b=coef(lm_summary)[1], L=L, adj.r.squared=lm_summary$adj.r.squared))
  return(list(lm_summary=lm_summary, used_data=used_data, a=a, b=b, L=L, adj.r.squared=lm_summary$adj.r.squared))
}

#'  calculate_stretched_exponential_model_by_optimization
#'  Formula to find linear fit between log(ln(L)-ln(s)) and log(N) that has the minimum value of 1-R**2 using optimization.
#'  @param results_dataframe Dataframe containing the data.
#'  @param y_parameter Parameter of the Y axis.
#'  @param x_parameter Parameter of the X axis.
#'  @param L_guess Initial values for the parameters to be optimized over.
#'  
calculate_stretched_exponential_model_by_optimization = function(results_dataframe, y_parameter, x_parameter, L_guess=c(2)){
  # Calculate mean of repetitions from same sample size
  results_dataframe_mean = results_dataframe %>% 
    arrange(get(x_parameter), rep) %>%
    group_by(get(x_parameter)) %>%
    summarise_at(vars(all_of(y_parameter)), list(mean=mean, median=median, sd=sd)) %>%
    rename(!!x_parameter := "get(x_parameter)") %>%
    ungroup()
  # Define s and N
  s = results_dataframe_mean$mean / max(results_dataframe_mean$mean)
  N = results_dataframe_mean$size
  # Define function to get linear fit
  linear_fit = function(data, par){
    # Calculate ln(L)-ln(s)
    L = par[1]
    s_rec = log(L) - log(data$s)
    # Calculate linear regression between ln(L)-ln(s) and ln(N)
    lm_result = summary(lm(log(s_rec)~log(data$N)))
    # The function returns the result of 1-R**2
    return(1-lm_result$adj.r.squared)
  }
  # Get minimization of error in linear fit
  #res = optim(par=L_guess, fn=linear_fit, data=data.frame(s=s, N=N), method="Nelder-Mead")
  res = optim(par=L_guess, fn=linear_fit, data=data.frame(s=s, N=N), method="Brent", lower=0, upper=5)
  L=res$par
  s_rec=log(L)-log(s)
  # This is the linear fit containing L
  lm_summary = summary(lm(log(s_rec)~log(N)))
  # Calculate alpha and b
  a = 1-coef(lm_summary)[2]
  b = exp((coef(lm_summary)[1] + log(a-1)))
  # This is the data used to obtain the final model
  used_data = data.frame(y=log(s_rec), x=log(N))
  return(list(lm_summary=lm_summary, used_data=used_data, a=a, b=b, L=L, adj.r.squared=lm_summary$adj.r.squared))
  #return(list(lm_summary=lm_summary, used_data=used_data, a=coef(lm_summary)[2], b=coef(lm_summary)[1], L=L, adj.r.squared=lm_summary$adj.r.squared))
}

#'  calculate_predictions_using_stretched_exponential_model_optimized
#'  Formula to calculate exponential decay of significant interactions from a list of sample sizes.
#'  This formula requires the use of the L parameter
#'  @param x List of sample sizes.
#'  @param a Slope coefficient.
#'  @param b Intercept coefficient.
#'  
calculate_predictions_using_stretched_exponential_model_optimized = function(x, L, a, b){
  y = L * exp((b*x**(-a+1))/(-a+1))
  #y = exp(log(L) - exp(b+a*log(x)))
  return(y)
}

#'  calculate_analytical_model
#'  Function to calculate the analytical model from different options
#'  @param results_dataframe Dataframe containing the data.
#'  @param y_parameter Parameter of the Y axis.
#'  @param x_parameter Parameter of the X axis.
#'  @param model Name of the model used.
#'  
calculate_analytical_model = function(results_dataframe, y_parameter, x_parameter, model, L=2){
  if(model == "Logarithmic"){
    model_output = get_logarithmic_tendency(results_dataframe=results_dataframe, y_parameter=y_parameter, x_parameter=x_parameter)
    model_result = (log(sort(unique(results_dataframe[[x_parameter]])))*model_output$a + model_output$b)
  } else if(model == "Stretched exponential (by optimization)"){
    model_output = calculate_stretched_exponential_model_by_optimization(results_dataframe=results_dataframe, y_parameter=y_parameter, x_parameter=x_parameter, L_guess=c(L))
    model_result = calculate_predictions_using_stretched_exponential_model_optimized(x=sort(unique(results_dataframe[[x_parameter]])), L=model_output$L, a=model_output$a, b=model_output$b)
  } else if(model == "Stretched exponential (by linear fit)"){
    model_output = calculate_stretched_exponential_model_by_linear_fit(results_dataframe=results_dataframe, y_parameter=y_parameter, x_parameter=x_parameter, L=L)
    model_result = calculate_predictions_using_stretched_exponential_model_optimized(x=sort(unique(results_dataframe[[x_parameter]])), L=model_output$L, a=model_output$a, b=model_output$b)
  } else if(model == "Stretched exponential (without L)"){
    model_output = calculate_stretched_exponential_model_without_L(results_dataframe=results_dataframe, y_parameter=y_parameter, x_parameter=x_parameter)
    model_result = calculate_predictions_using_stretched_exponential_model_without_L(x=sort(unique(results_dataframe[[x_parameter]])), a=model_output$a, b=model_output$b)
  }
  return(list(model_output=model_output, model_result=model_result))
}

#'  calculate_prediction_from_analytical_model
#'  Function to calculate predictions using a specific analytical model
#'  @param model Name of the model used.
#'  @param x_list List of values in the x axis (e.g. size)
#'  @param a Parameter a.
#'  @param b Parameter b.
#'  @param L Parameter L.
#'  
calculate_prediction_from_analytical_model = function(model, x_list, a, b, L){
  if(model == "Logarithmic"){
    prediction_result = (log(x_list)*a + b)
  } else if((model == "Stretched exponential (by optimization)") | (model == "Stretched exponential (by linear fit)")){
    prediction_result = calculate_predictions_using_stretched_exponential_model_optimized(x=x_list, L=L, a=a, b=b)
  } else if(model == "Stretched exponential (without L)"){
    prediction_result = calculate_predictions_using_stretched_exponential_model_without_L(x=x_list, a=a, b=b)
  }
  return(prediction_result)
}

#'  get_formula
#'  Function get formula from parameters
#'  @param model Name of the model used.
#'  @param a Parameter a.
#'  @param b Parameter b.
#'  @param L Parameter L.
#'  
get_formula = function(model, a, b, L){
  if(model == "Logarithmic"){
    formula_name = paste("F(x) = (", formatC(round(a, 2), format = "e", digits = 2), ")*ln(x) + (", formatC(round(b, 2), format = "e", digits = 2), ")", sep="")
  } else if((model == "Stretched exponential (by optimization)") | (model == "Stretched exponential (by linear fit)")){
    formula_name = paste("F(x) = ", round(L, 2), " * exp[(", round(b, 2), " * x**((", round(a, 2), ") + 1)) / ((", round(a, 2), ") + 1) ]", sep="")
  } else if(model == "Stretched exponential (without L)"){
    formula_name = paste("F(x) = exp[ (exp(", round(b, 2), ") * x**(1 + (", round(a, 2), "))) / (1 + (", round(a, 2), ")) ]", sep="")
  }
  return(formula_name)
}


#-----------#
# Read data #
#-----------#

topology_results_df = fread(topology_results_file)
ppi_results_df = fread(ppi_results_file)
disease_genes_results_df = fread(disease_genes_results_file)
essential_genes_results_df = fread(essential_genes_results_file) %>% rename("num_essential_components" = "num_components", "num_essential_lcc_nodes" = "num_lcc_nodes", "num_essential_lcc_edges" = "num_lcc_edges", "essential_lcc_z" = "lcc_z", "essential_lcc_pvalue" = "lcc_pvalue")
results_df = inner_join(topology_results_df, ppi_results_df, by = c("method", "dataset", "type_dataset", "size", "rep", "type_correlation", "threshold")) %>% inner_join(disease_genes_results_df, by = c("method", "dataset", "type_dataset", "size", "rep", "type_correlation", "threshold")) %>% inner_join(essential_genes_results_df, by = c("method", "dataset", "type_dataset", "size", "rep", "type_correlation", "threshold"))
results_df$type_dataset = paste(results_df$dataset, results_df$type_dataset, sep=":") # Join dataset and type_dataset
results_df$type_dataset = tolower(results_df$type_dataset)
results_df$threshold = as.character(results_df$threshold) # Consider threshold as a discrete variable
results_df$log_disease_lcc_pvalue = abs(log10(results_df$disease_lcc_pvalue))
results_df$log_disease_lcc_pvalue = replace(results_df$log_disease_lcc_pvalue, is.infinite(results_df$log_disease_lcc_pvalue),abs(log(0.000001))) # Replace infinite values by very high values
numbers_complete_graph_df = fread(numbers_complete_graph_file) %>% rename("type_dataset" = "dataset")
numbers_complete_graph_df$type_dataset = tolower(numbers_complete_graph_df$type_dataset)

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

  observeEvent(input$topology_analytical, {
    
    freezeReactiveValue(input, "topology_type_normalization") # Freeze topology_type_normalization until updateSelectInput is completed
    if ((isTRUE(input$topology_analytical)) & (("Stretched exponential (by optimization)" %in% input$topology_type_analytical_model) | ("Stretched exponential (by linear fit)" %in% input$topology_type_analytical_model))){
      updateSelectInput(session = session, 
                        inputId = "topology_type_normalization", 
                        label = "Type of y axis normalization:", 
                        choices = c("Divide by max. value" = "divide.max.value", "Divide by max. possible value" = "divide.max.possible.value", "Divide by value of convergence" = "divide.L"),
                        selected = "divide.max.value")
    } else {
      updateSelectInput(session = session, 
                        inputId = "topology_type_normalization", 
                        label = "Type of y axis normalization:", 
                        choices = c("Divide by max. value" = "divide.max.value", "Divide by max. possible value" = "divide.max.possible.value"),
                        selected = "divide.max.value")
    }
  })
  
  
  process_inputs = reactive({
  #output$topologyBoxPlot <- renderPlot({
    
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
      type_datasets_gtex = paste("gtex", type_datasets_gtex, sep=":")
    }
    type_datasets_scipher = c()
    if ("scipher" %in% c(input$topology_dataset)){
      type_datasets_scipher=c(input$topology_type_dataset_scipher)
      type_datasets_scipher = paste("scipher", type_datasets_scipher, sep=":")
    }
    type_datasets_tcga = c()
    if ("tcga" %in% c(input$topology_dataset)){
      type_datasets_tcga=c(input$topology_tcga_project)
      type_datasets_tcga = paste("tcga", type_datasets_tcga, sep=":")
    }
    type_datasets = c(type_datasets_gtex, type_datasets_scipher, type_datasets_tcga)
      
    #--------------------#
    # Select information #
    #--------------------#

    results_selected_df = results_df %>% filter((type_dataset %in% type_datasets) & (method %in% c(input$topology_method)) & (type_correlation %in% c(input$topology_type_correlation)) & (threshold %in% c(input$topology_pvalue_threshold)))
    #results_selected_df = results_df %>% filter((type_dataset %in% c("tcga:tcga")) & (method %in% c("pearson")) & (type_correlation %in% c("all")) & (threshold %in% c(0.05))) # For a test
    #results_selected_df = results_df %>% filter((type_dataset %in% c("tcga:tcga-brca", "tcga:tcga-ucec")) & (method %in% c("pearson")) & (type_correlation %in% c("all")) & (threshold %in% c(0.05))) # For a test with multiple parameters

    # Merge with total number of edges
    results_selected_df = results_selected_df %>% inner_join(numbers_complete_graph_df %>% select("type_dataset", "total_num_edges"), by = "type_dataset") %>% group_by_at(c("dataset", "type_dataset")) %>% mutate(total_num_edges_norm=total_num_edges/max(num_edges)) %>% ungroup() %>% as.data.frame()
    
    # Select by diseases
    group_vars = c("type_dataset", "method", "size", "type_correlation", "threshold", "total_num_edges")
    if(input$type_analysis == "disease_genes"){
      diseases_selected <- unlist(strsplit(renderText(input$diseases, sep='___')(), split='___'))
      results_selected_df %<>% filter(disease %in% diseases_selected) %>% select(-disease_class) %>% unique()
      group_vars = c(group_vars, "disease")
    } else {
      results_selected_df %<>% select(!(c("num_disease_genes", "num_disease_edges", "fraction_disease_genes", "num_disease_components", "disease", "disease_class", "num_disease_lcc_nodes", "num_disease_lcc_edges", "fraction_disease_lcc_nodes", "disease_lcc_z", "disease_lcc_pvalue", "log_disease_lcc_pvalue"))) %>% unique()
    }
    
    
    #--------------------------#
    # Calculate mean / sd line #
    #--------------------------#
    
    # Calculate mean and standard deviation tendency
    topology_results_selected_by_size_df = results_selected_df %>%
      group_by_at(group_vars) %>%
      summarise_at(vars(all_of(input$boxplot_parameter)), list(mean=mean, sd=sd)) %>%
      arrange(size)
    
    # Check if the user wants to plot mean or standard deviation
    if(isTRUE(input$topology_sd)){
      combination_metric = "sd"
    } else{
      combination_metric = "mean"
    }
    
    # Check which parameter should be plotted multiple times (if there are multiple checkboxes with multiple elements selected)
    # In this case, the priority is given by the following list (type_dataset > method...)
    multiple_options_params = c("type_dataset", "method", "type_correlation", "threshold", "disease")
    selected_fill_parameter = NULL
    for (multiple_options_param in multiple_options_params){
      if (length(unique(results_selected_df[[multiple_options_param]])) > 1){
        selected_fill_parameter = multiple_options_param
        break
      }
    }

    #----------------------------#
    # Calculate analytical curve #
    #----------------------------#
    
    # Check if user wants to plot analytical model
    cols = c("model", "model_result", "size", "a", "b", "L", "adj.r.squared", "max_value_in_dataset", "formula", "type_dataset", "fill_parameter")
    topology_results_selected_analytical_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
    cols = c("model", "y", "x", "regression", "type_dataset", "fill_parameter")
    stretched_exponential_regression_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
    cols = c("model", "model_result", "size", "max_value_in_dataset", "type_dataset", "fill_parameter")
    predicted_results_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
    analytical_model_summary_df = data.frame()
    N_vals = seq(10, 50000, 10)
    
    types_analytical_model = c("Stretched exponential (by optimization)", "Stretched exponential (by linear fit)", "Stretched exponential (without L)", "Logarithmic")
    if (isTruthy(input$topology_analytical)){
      for(model in types_analytical_model){
        if (is.null(selected_fill_parameter)){
          # Calculate the analytical model
          output_variable = calculate_analytical_model(results_dataframe=results_selected_df, y_parameter=input$boxplot_parameter, x_parameter="size", model=model, L=unique(results_selected_df$total_num_edges_norm))
          # Calculate prediction
          prediction_result = calculate_prediction_from_analytical_model(model=model, x_list=N_vals, a=output_variable$model_output$a, b=output_variable$model_output$b, L=output_variable$model_output$L)
          # Un-normalize y axis of analytical curve (if necessary), because the analytical curve of stretched exponential is normalized by default
          max_value_in_dataset = max(results_selected_df[[input$boxplot_parameter]])
          if (!(model == "Logarithmic")){
            output_variable$model_result = output_variable$model_result * max_value_in_dataset
            prediction_result = prediction_result * max_value_in_dataset
            # Save data used to create the model
            stretched_exponential_regression_df = rbind(stretched_exponential_regression_df, cbind(data.frame(model=model), output_variable$model_output$used_data, data.frame(regression_line=(output_variable$model_output$a * output_variable$model_output$used_data$x + output_variable$model_output$b), type_dataset=unique(results_selected_df$type_dataset)[1])))
          }
          # Calculate formula
          formula = get_formula(model=model, a=output_variable$model_output$a, b=output_variable$model_output$b, L=output_variable$model_output$L)
          # Add results
          topology_results_selected_analytical_df = rbind(topology_results_selected_analytical_df, data.frame(model=model, model_result=output_variable$model_result, size=sort(unique(results_selected_df$size)), a=output_variable$model_output$a, b=output_variable$model_output$b, L=output_variable$model_output$L, adj.r.squared=output_variable$model_output$lm_summary$adj.r.squared, max_value_in_dataset=max_value_in_dataset, formula=formula, type_dataset=unique(results_selected_df$type_dataset)[1]))
          predicted_results_df = rbind(predicted_results_df, data.frame(model=model, model_result=prediction_result, size=N_vals, max_value_in_dataset=max_value_in_dataset, type_dataset=unique(results_selected_df$type_dataset)[1]))
        } else{
          for(selected_parameter in unique(results_selected_df[[selected_fill_parameter]])){
            # Select results for a specific parameter
            results_selected_by_parameter_df = results_selected_df %>% filter(!!as.symbol(selected_fill_parameter) == selected_parameter)
            # Calculate the analytical model
            output_variable = calculate_analytical_model(results_dataframe=results_selected_by_parameter_df, y_parameter=input$boxplot_parameter, x_parameter="size", model=model, L=unique(results_selected_by_parameter_df$total_num_edges_norm))
            # Calculate prediction
            prediction_result = calculate_prediction_from_analytical_model(model=model, x_list=N_vals, a=output_variable$model_output$a, b=output_variable$model_output$b, L=output_variable$model_output$L)
            # Un-normalize y axis of analytical curve (if necessary), because the analytical curve of stretched exponential is normalized by default
            max_value_in_dataset = max(results_selected_by_parameter_df[[input$boxplot_parameter]])
            if (!(model == "Logarithmic")){
              output_variable$model_result = output_variable$model_result * max_value_in_dataset
              prediction_result = prediction_result * max_value_in_dataset
              # Save data used to create the model
              stretched_exponential_regression_df = rbind(stretched_exponential_regression_df, cbind(data.frame(model=model), output_variable$model_output$used_data, data.frame(regression_line=(output_variable$model_output$a * output_variable$model_output$used_data$x + output_variable$model_output$b), type_dataset=unique(results_selected_by_parameter_df$type_dataset)[1], fill_parameter=selected_parameter)))
            }
            # Calculate formula
            formula = get_formula(model=model, a=output_variable$model_output$a, b=output_variable$model_output$b, L=output_variable$model_output$L)
            # Add results
            topology_results_selected_analytical_df = rbind(topology_results_selected_analytical_df, data.frame(model=model, model_result=output_variable$model_result, size=sort(unique(results_selected_by_parameter_df$size)), a=output_variable$model_output$a, b=output_variable$model_output$b, L=output_variable$model_output$L, adj.r.squared=output_variable$model_output$lm_summary$adj.r.squared, max_value_in_dataset=max_value_in_dataset, formula=formula, type_dataset=unique(results_selected_by_parameter_df$type_dataset), fill_parameter=selected_parameter))
            predicted_results_df = rbind(predicted_results_df, data.frame(model=model, model_result=prediction_result, size=N_vals, max_value_in_dataset=max_value_in_dataset, type_dataset=unique(results_selected_by_parameter_df$type_dataset), fill_parameter=selected_parameter))
          }
        }
      }
      # Rename fill parameter if necessary
      if(("fill_parameter" %in% colnames(predicted_results_df)) && (!(is.null(selected_fill_parameter))) && (!(selected_fill_parameter == "type_dataset"))){
        topology_results_selected_analytical_df = topology_results_selected_analytical_df %>% rename(!!selected_fill_parameter := fill_parameter)
        predicted_results_df = predicted_results_df %>% rename(!!selected_fill_parameter := fill_parameter)
        stretched_exponential_regression_df = stretched_exponential_regression_df %>% rename(!!selected_fill_parameter := fill_parameter)
      }
      # Join predictions with real data to calculate error
      topology_results_selected_pred_vs_reality_df = topology_results_selected_analytical_df %>% 
        inner_join(topology_results_selected_by_size_df %>% ungroup() %>% select("size", !!combination_metric, !!selected_fill_parameter), by=c("size", selected_fill_parameter)) %>% 
        unique()
      # Calculate relative error
      topology_results_selected_pred_vs_reality_df$relative.error = abs((topology_results_selected_pred_vs_reality_df[[combination_metric]] - topology_results_selected_pred_vs_reality_df$model_result)) / topology_results_selected_pred_vs_reality_df[[combination_metric]]
      # Create summary table
      analytical_model_summary_df = topology_results_selected_pred_vs_reality_df %>% 
        filter(model %in% c(input$topology_type_analytical_model)) %>% 
        select(one_of(selected_fill_parameter), "type_dataset", "model", "a", "b", "L", "max_value_in_dataset", "formula", "adj.r.squared", "relative.error") %>% 
        unique() %>%
        # Calculate relative error mean
        group_by_at(c(selected_fill_parameter, "model", "formula", "adj.r.squared")) %>%
        mutate(relative.error.mean=mean(relative.error)) %>%
        select(!("relative.error")) %>%
        unique() %>%
        ungroup()
      analytical_model_summary_df$unnorm_L = analytical_model_summary_df$L * analytical_model_summary_df$max_value_in_dataset
      analytical_model_summary_df$a_b_ratio = analytical_model_summary_df$a / analytical_model_summary_df$b
      analytical_model_summary_df = analytical_model_summary_df %>% inner_join(numbers_complete_graph_df %>% select("type_dataset", "total_num_edges"), by = "type_dataset") %>% group_by_at(selected_fill_parameter) %>% mutate(density = unnorm_L/total_num_edges)
      # Render table
      #output$topologyAnalyticalModelTable = renderTable(analytical_model_summary_df)
      output$topologyAnalyticalModelTable = DT::renderDataTable(analytical_model_summary_df)
    } else {
      # Rename fill parameter if necessary
      if(("fill_parameter" %in% colnames(predicted_results_df)) && (!(is.null(selected_fill_parameter))) && (!(selected_fill_parameter == "type_dataset"))){
        topology_results_selected_analytical_df = topology_results_selected_analytical_df %>% rename(!!selected_fill_parameter := fill_parameter)
      }
      # Remove table if it is rendered
      removeUI(selector = "#topologyAnalyticalModelTableID")
    }

    #------------------------------------#
    # Normalize data if required by user #
    #------------------------------------#
    
    # Re-scale X axis
    if (isTRUE(input$topology_normalize_x)){
      results_selected_df = results_selected_df %>% group_by_at(selected_fill_parameter) %>% mutate(max_size = max(size)) %>% mutate(norm = size/max_size) %>% select(!(all_of(c("size", "max_size")))) %>% rename(size = norm)
      topology_results_selected_by_size_df = topology_results_selected_by_size_df %>% group_by_at(selected_fill_parameter) %>% mutate(max_size=max(size)) %>% mutate(norm = size/max_size) %>% select(!(all_of(c("size", "max_size")))) %>% rename(size = norm)
      topology_results_selected_analytical_df = topology_results_selected_analytical_df %>% group_by_at(selected_fill_parameter) %>% mutate(max_size = max(size)) %>% mutate(norm = size/max_size) %>% select(!(all_of(c("size", "max_size")))) %>% rename(size = norm)
      predicted_results_df = predicted_results_df %>% group_by_at(selected_fill_parameter) %>% mutate(max_size = max(size)) %>% mutate(norm = size/max_size) %>% select(!(all_of(c("size", "max_size")))) %>% rename(size = norm)
    }
    
    # Rescale Y axis
    if(isTRUE(input$topology_normalize_y)){
      if(input$topology_type_normalization == "divide.max.value"){
        results_selected_df = results_selected_df %>% group_by_at(selected_fill_parameter) %>% mutate(max_parameter = max(get(input$boxplot_parameter))) %>% mutate(norm = get(input$boxplot_parameter)/max_parameter) %>% rename(unnorm = input$boxplot_parameter) %>% rename(!!input$boxplot_parameter := norm)
        topology_results_selected_by_size_df = topology_results_selected_by_size_df %>% group_by_at(selected_fill_parameter) %>% mutate(max_mean=max(mean)) %>% mutate(mean_norm=mean/max_mean) %>% mutate(max_sd=max(sd)) %>% mutate(sd_norm=sd/max_sd)
        topology_results_selected_analytical_df = topology_results_selected_analytical_df %>% group_by_at(selected_fill_parameter) %>% mutate(model_result = model_result/max_value_in_dataset)
        predicted_results_df = predicted_results_df %>% group_by_at(selected_fill_parameter) %>% mutate(model_result = model_result/max_value_in_dataset)
      } else if(input$topology_type_normalization == "divide.L"){
        if (is.null(selected_fill_parameter)){
          results_selected_df = results_selected_df %>% inner_join((topology_results_selected_analytical_df %>% filter(model %in% input$topology_type_analytical_model) %>% select("L", "max_value_in_dataset", "type_dataset") %>% unique()), by="type_dataset") %>% mutate(norm = (get(input$boxplot_parameter)/max_value_in_dataset)/L) %>% rename(unnorm = input$boxplot_parameter) %>% rename(!!input$boxplot_parameter := norm)
          topology_results_selected_by_size_df = topology_results_selected_by_size_df %>% inner_join((topology_results_selected_analytical_df %>% filter(model %in% input$topology_type_analytical_model) %>% select("L", "max_value_in_dataset", "type_dataset") %>% unique()), by="type_dataset") %>% group_by_at(selected_fill_parameter) %>% mutate(mean_norm = (mean/max_value_in_dataset)/L)
          predicted_results_df = predicted_results_df %>% inner_join((topology_results_selected_analytical_df %>% filter(model %in% input$topology_type_analytical_model) %>% select("L", "type_dataset") %>% unique()), by="type_dataset") %>% group_by_at(selected_fill_parameter) %>% mutate(model_result = (model_result / max_value_in_dataset)/L)
          topology_results_selected_analytical_df = topology_results_selected_analytical_df %>% group_by_at(selected_fill_parameter) %>% mutate(norm = (model_result/max_value_in_dataset)/L) %>% rename(unnorm = model_result) %>% rename(model_result = norm)
        } else{
          results_selected_df = results_selected_df %>% inner_join((topology_results_selected_analytical_df %>% filter(model %in% input$topology_type_analytical_model) %>% select("L", "max_value_in_dataset", !!selected_fill_parameter) %>% unique()), by=selected_fill_parameter) %>% group_by_at(selected_fill_parameter) %>% mutate(norm = (get(input$boxplot_parameter)/max_value_in_dataset)/L) %>% rename(unnorm = input$boxplot_parameter) %>% rename(!!input$boxplot_parameter := norm)
          topology_results_selected_by_size_df = topology_results_selected_by_size_df %>% inner_join((topology_results_selected_analytical_df %>% filter(model %in% input$topology_type_analytical_model) %>% select("L", "max_value_in_dataset", !!selected_fill_parameter) %>% unique()), by=selected_fill_parameter) %>% group_by_at(selected_fill_parameter) %>% mutate(mean_norm = (mean/max_value_in_dataset)/L)
          predicted_results_df = predicted_results_df %>% inner_join((topology_results_selected_analytical_df %>% filter(model %in% input$topology_type_analytical_model) %>% select("L", !!selected_fill_parameter) %>% unique()), by=selected_fill_parameter) %>% group_by_at(selected_fill_parameter) %>% mutate(model_result = (model_result / max_value_in_dataset)/L)
          topology_results_selected_analytical_df = topology_results_selected_analytical_df %>% group_by_at(selected_fill_parameter) %>% mutate(norm = (model_result/max_value_in_dataset)/L) %>% rename(unnorm = model_result) %>% rename(model_result = norm)
        }
      } else if(input$topology_type_normalization == "divide.max.possible.value"){
        results_selected_df = results_selected_df %>% 
          #inner_join(numbers_complete_graph_df, by = "type_dataset") %>% 
          group_by_at(selected_fill_parameter) %>% mutate(norm = get(input$boxplot_parameter)/total_num_edges) %>% rename(unnorm = input$boxplot_parameter) %>% rename(!!input$boxplot_parameter := norm)
        topology_results_selected_by_size_df = topology_results_selected_by_size_df %>% inner_join(numbers_complete_graph_df, by = "type_dataset") %>% group_by_at(selected_fill_parameter) %>% mutate(mean_norm=mean/total_num_edges)
        if (isTruthy(input$topology_analytical)){
          topology_results_selected_analytical_df = topology_results_selected_analytical_df %>% inner_join(numbers_complete_graph_df, by = "type_dataset") %>% group_by_at(selected_fill_parameter) %>% mutate(norm = model_result/total_num_edges) %>% rename(unnorm = model_result) %>% rename(model_result = norm)
          predicted_results_df = predicted_results_df %>% inner_join(numbers_complete_graph_df, by = "type_dataset") %>% group_by_at(selected_fill_parameter) %>% mutate(norm = model_result/total_num_edges) %>% rename(unnorm = model_result) %>% rename(model_result = norm)
          # print size for percentage of edges
          for(type_dataset_selected in unique(predicted_results_df$type_dataset)){
            predicted_results_specific = predicted_results_df %>% filter((model %in% c(input$topology_type_analytical_model)) & (type_dataset==type_dataset_selected))
            #print(predicted_results_specific[which.min(abs(0.25-predicted_results_specific$model_result)),])
            #print(predicted_results_specific[which.min(abs(0.50-predicted_results_specific$model_result)),])
            #print(predicted_results_specific[which.min(abs(0.75-predicted_results_specific$model_result)),])
          }
        }
      }
      combination_metric = paste(combination_metric, "norm", sep="_")
    }

    # if (isTruthy(input$topology_analytical)){
    #   if (is.null(selected_fill_parameter)){
    #     # Calculate logarithmic fit for the whole selected data
    #     log_results = get_logarithmic_tendency(results_dataframe=results_selected_df, y_parameter=input$boxplot_parameter, x_parameter="size")
    #     dec_model = calculate_stretched_exponential_model_without_L(results_dataframe=results_selected_df, y_parameter=input$boxplot_parameter, x_parameter="size")
    #     dec_model_results = calculate_predictions_using_stretched_exponential_model_without_L(x=sort(unique(results_selected_df$size)), a=coef(dec_model$lm_summary)[2], b=coef(dec_model$lm_summary)[1])
    #     dec_smax_model = calculate_stretched_exponential_model_by_optimization(results_dataframe=results_selected_df, y_parameter=input$boxplot_parameter, x_parameter="size", L_guess=c(2))
    #     dec_smax_model_results = calculate_predictions_using_stretched_exponential_model_optimized(x=sort(unique(results_selected_df$size)), L=dec_smax_model$L, a=coef(dec_smax_model$lm_summary)[2], b=coef(dec_smax_model$lm_summary)[1])
    #     if(input$topology_type_analytical_model == "Stretched exponential (by optimization)"){
    #       stretched_exponential_regression_df = cbind(data.frame(model="Stretched exponential (by optimization)"), dec_smax_model$dat, data.frame(regression_line=(coef(dec_smax_model$lm_summary)[2] * dec_smax_model$dat$x + coef(dec_smax_model$lm_summary)[1])))
    #     }
    #     if(input$topology_type_analytical_model == "Stretched exponential (without L)"){
    #       stretched_exponential_regression_df = cbind(data.frame(model="Stretched exponential (without L)"), dec_model$dat, data.frame(regression_line=(coef(dec_model$lm_summary)[2] * dec_model$dat$x + coef(dec_model$lm_summary)[1])))
    #     }
    #     topology_results_selected_analytical_df = rbind(topology_results_selected_analytical_df, data.frame(model="Logarithmic", model_result=(log(sort(unique(results_selected_df$size)))*coef(log_results$lm_summary)[2] + coef(log_results$lm_summary)[1]), size=sort(unique(results_selected_df$size)), formula=log_results$formula_name, adj.r.squared=log_results$lm_summary$adj.r.squared))
    #     topology_results_selected_analytical_df = rbind(topology_results_selected_analytical_df, data.frame(model="Stretched exponential (without L)", model_result=dec_model_results, size=sort(unique(results_selected_df$size)), formula=dec_model$formula_name, adj.r.squared=dec_model$lm_summary$adj.r.squared))
    #     topology_results_selected_analytical_df = rbind(topology_results_selected_analytical_df, data.frame(model="Stretched exponential (by optimization)", model_result=dec_smax_model_results, size=sort(unique(results_selected_df$size)), formula=dec_smax_model$formula_name, adj.r.squared=dec_smax_model$lm_summary$adj.r.squared))
    #     # Calculate predictions
    #     dec_model_predictions = calculate_predictions_using_stretched_exponential_model_without_L(x=N_vals, a=coef(dec_model$lm_summary)[2], b=coef(dec_model$lm_summary)[1])
    #     dec_smax_model_predictions = calculate_predictions_using_stretched_exponential_model_optimized(x=N_vals, L=dec_smax_model$L, a=coef(dec_smax_model$lm_summary)[2], b=coef(dec_smax_model$lm_summary)[1])
    #     predicted_results_df = rbind(predicted_results_df, data.frame(model="Logarithmic", model_result=(N_vals*coef(log_results$lm_summary)[2] + coef(log_results$lm_summary)[1]), size=N_vals, type_dataset=unique(results_selected_df$type_dataset)[1]))
    #     predicted_results_df = rbind(predicted_results_df, data.frame(model="Stretched exponential (without L)", model_result=dec_model_predictions, size=N_vals, type_dataset=unique(results_selected_df$type_dataset)[1]))
    #     predicted_results_df = rbind(predicted_results_df, data.frame(model="Stretched exponential (by optimization)", model_result=dec_smax_model_predictions, size=N_vals, type_dataset=unique(results_selected_df$type_dataset)[1]))
    #     # Un-normalize y axis of analytical curve if necessary (because the analyticsal curve of Stretched exponential (by optimization) is normalized by default)
    #     if(!(isTRUE(input$topology_normalize_y))){
    #       topology_results_selected_analytical_df$model_result = ifelse(topology_results_selected_analytical_df$model %in% c("Stretched exponential (without L)", "Stretched exponential (by optimization)"), topology_results_selected_analytical_df$model_result * max(results_selected_df[[input$boxplot_parameter]]), topology_results_selected_analytical_df$model_result)
    #       predicted_results_df$model_result = ifelse(predicted_results_df$model %in% c("Stretched exponential (without L)", "Stretched exponential (by optimization)"), predicted_results_df$model_result * max(results_selected_df[[input$boxplot_parameter]]), predicted_results_df$model_result)
    #     }
    #     # Un-normalize and re-normalize by total num of edges in complete graph
    #     if((isTRUE(input$topology_normalize_y)) & (input$topology_type_normalization == "divide.max.possible.value")){
    #       topology_results_selected_analytical_df$model_result = ifelse(topology_results_selected_analytical_df$model %in% c("Stretched exponential (without L)", "Stretched exponential (by optimization)"), topology_results_selected_analytical_df$model_result * max(results_selected_df$unnorm) / unique(results_selected_df$total_num_edges), topology_results_selected_analytical_df$model_result)
    #       predicted_results_df = predicted_results_df %>% inner_join((results_selected_df %>% select("type_dataset", "total_num_edges") %>% unique()), by=c("type_dataset")) %>% mutate(norm = model_result * max(results_selected_df$unnorm)/total_num_edges) %>% select(!c("model_result", "total_num_edges")) %>% rename(model_result = norm)
    #     }
    #     # Divide by L to re-scale the axis using as maximum the L value
    #     if((isTRUE(input$topology_normalize_y)) & (input$topology_type_normalization == "divide.L")){
    #       topology_results_selected_analytical_df$model_result = topology_results_selected_analytical_df$model_result / dec_smax_model$L
    #       predicted_results_df$model_result = predicted_results_df$model_result / dec_smax_model$L
    #       results_selected_df = results_selected_df %>% mutate(norm = get(input$boxplot_parameter)/dec_smax_model$L) %>% select(!(all_of(c(!!input$boxplot_parameter)))) %>% rename(!!input$boxplot_parameter := norm)
    #       topology_results_selected_by_size_df$mean_norm = topology_results_selected_by_size_df$mean_norm / dec_smax_model$L
    #     }
    #   } else {
    #     # Calculate logarithmic fit for each selected fill parameter (e.g. each dataset)
    #     for (selected_parameter in unique(results_selected_df[[selected_fill_parameter]])){
    #       results_selected_by_parameter_df = results_selected_df %>% filter(!!as.symbol(selected_fill_parameter) == selected_parameter)
    #       log_results = get_logarithmic_tendency(results_dataframe=results_selected_by_parameter_df, y_parameter=input$boxplot_parameter, x_parameter="size")
    #       dec_model = calculate_stretched_exponential_model_without_L(results_dataframe=results_selected_by_parameter_df, y_parameter=input$boxplot_parameter, x_parameter="size")
    #       dec_model_results = calculate_predictions_using_stretched_exponential_model_without_L(x=sort(unique(results_selected_by_parameter_df$size)), a=coef(dec_model$lm_summary)[2], b=coef(dec_model$lm_summary)[1])
    #       dec_smax_model = calculate_stretched_exponential_model_by_optimization(results_dataframe=results_selected_by_parameter_df, y_parameter=input$boxplot_parameter, x_parameter="size", L_guess=c(2))
    #       dec_smax_model_results = calculate_predictions_using_stretched_exponential_model_optimized(x=sort(unique(results_selected_by_parameter_df$size)), L=dec_smax_model$L, a=coef(dec_smax_model$lm_summary)[2], b=coef(dec_smax_model$lm_summary)[1])
    #       if(input$topology_type_analytical_model == "Stretched exponential (by optimization)"){
    #         stretched_exponential_regression_df = rbind(stretched_exponential_regression_df, cbind(data.frame(model="Stretched exponential (by optimization)"), dec_smax_model$dat, data.frame(regression_line=(coef(dec_smax_model$lm_summary)[2] * dec_smax_model$dat$x + coef(dec_smax_model$lm_summary)[1])), (data.frame(parameter = selected_parameter) %>% rename(!!selected_fill_parameter := parameter))))
    #       }
    #       if(input$topology_type_analytical_model == "Stretched exponential (without L)"){
    #         stretched_exponential_regression_df = rbind(stretched_exponential_regression_df, cbind(data.frame(model="Stretched exponential (without L)"), dec_model$dat, data.frame(regression_line=(coef(dec_model$lm_summary)[2] * dec_model$dat$x + coef(dec_model$lm_summary)[1])), (data.frame(parameter = selected_parameter) %>% rename(!!selected_fill_parameter := parameter))))
    #       }
    #       topology_results_selected_analytical_df = rbind(topology_results_selected_analytical_df, data.frame(model="Logarithmic", model_result=(log(sort(unique(results_selected_by_parameter_df$size)))*coef(log_results$lm_summary)[2] + coef(log_results$lm_summary)[1]), size=sort(unique(results_selected_by_parameter_df$size)), formula=log_results$formula_name, adj.r.squared=log_results$lm_summary$adj.r.squared, parameter=selected_parameter) %>% rename(!!selected_fill_parameter := parameter))
    #       topology_results_selected_analytical_df = rbind(topology_results_selected_analytical_df, data.frame(model="Stretched exponential (without L)", model_result=dec_model_results, size=sort(unique(results_selected_by_parameter_df$size)), formula=dec_model$formula_name, adj.r.squared=dec_model$lm_summary$adj.r.squared, parameter=selected_parameter) %>% rename(!!selected_fill_parameter := parameter))
    #       topology_results_selected_analytical_df = rbind(topology_results_selected_analytical_df, data.frame(model="Stretched exponential (by optimization)", model_result=dec_smax_model_results, size=sort(unique(results_selected_by_parameter_df$size)), formula=dec_smax_model$formula_name, adj.r.squared=dec_smax_model$lm_summary$adj.r.squared, parameter=selected_parameter) %>% rename(!!selected_fill_parameter := parameter))
    #       # Calculate predictions
    #       dec_model_predictions = calculate_predictions_using_stretched_exponential_model_without_L(x=N_vals, a=coef(dec_model$lm_summary)[2], b=coef(dec_model$lm_summary)[1])
    #       dec_smax_model_predictions = calculate_predictions_using_stretched_exponential_model_optimized(x=N_vals, L=dec_smax_model$L, a=coef(dec_smax_model$lm_summary)[2], b=coef(dec_smax_model$lm_summary)[1])
    #       predicted_results_df = rbind(predicted_results_df, data.frame(model="Logarithmic", model_result=(log(N_vals)*coef(log_results$lm_summary)[2] + coef(log_results$lm_summary)[1]), size=N_vals, type_dataset=unique(results_selected_by_parameter_df$type_dataset), fill_parameter=selected_parameter))
    #       predicted_results_df = rbind(predicted_results_df, data.frame(model="Stretched exponential (without L)", model_result=dec_model_predictions, size=N_vals, type_dataset=unique(results_selected_by_parameter_df$type_dataset), fill_parameter=selected_parameter))
    #       predicted_results_df = rbind(predicted_results_df, data.frame(model="Stretched exponential (by optimization)", model_result=dec_smax_model_predictions, size=N_vals, type_dataset=unique(results_selected_by_parameter_df$type_dataset), fill_parameter=selected_parameter))
    #       # Un-normalize y axis of analytical curve if necessary (because the analytical curve of Stretched exponential (by optimization) is normalized by default)
    #       if(!(isTRUE(input$topology_normalize_y))){
    #         topology_results_selected_analytical_df$model_result = ifelse((topology_results_selected_analytical_df[[selected_fill_parameter]] == selected_parameter & topology_results_selected_analytical_df$model %in% c("Stretched exponential (without L)", "Stretched exponential (by optimization)")), topology_results_selected_analytical_df$model_result * max(results_selected_by_parameter_df[[input$boxplot_parameter]]), topology_results_selected_analytical_df$model_result)
    #         predicted_results_df$model_result = ifelse((predicted_results_df$fill_parameter == selected_parameter & predicted_results_df$model %in% c("Stretched exponential (without L)", "Stretched exponential (by optimization)")), predicted_results_df$model_result * max(results_selected_by_parameter_df[[input$boxplot_parameter]]), predicted_results_df$model_result)
    #       }
    #       # Un-normalize and re-normalize by total num of edges in complete graph
    #       if((isTRUE(input$topology_normalize_y)) & (input$topology_type_normalization == "divide.max.possible.value")){
    #         topology_results_selected_analytical_df$model_result = ifelse((topology_results_selected_analytical_df[[selected_fill_parameter]] == selected_parameter & topology_results_selected_analytical_df$model %in% c("Stretched exponential (without L)", "Stretched exponential (by optimization)")), topology_results_selected_analytical_df$model_result * max(results_selected_by_parameter_df$unnorm) / unique(results_selected_by_parameter_df$total_num_edges), topology_results_selected_analytical_df$model_result)
    #         # Join predicted results with total number of edges
    #         predicted_results_df = predicted_results_df %>% inner_join((results_selected_df %>% select("type_dataset", "total_num_edges") %>% unique()), by=c("type_dataset"))
    #         # Normalize result for selected parameter
    #         predicted_results_df$model_result = ifelse((predicted_results_df$fill_parameter == selected_parameter & predicted_results_df$model %in% c("Stretched exponential (without L)", "Stretched exponential (by optimization)")), predicted_results_df$model_result * max(results_selected_by_parameter_df$unnorm) / predicted_results_df$total_num_edges, predicted_results_df$model_result)
    #         # Remove total number of edges
    #         predicted_results_df = predicted_results_df %>% select(!c("total_num_edges"))
    #       }
    #       # Divide by L to re-scale the axis using as maximum the L value
    #       if((isTRUE(input$topology_normalize_y)) & (input$topology_type_normalization == "divide.L")){
    #         topology_results_selected_analytical_df$model_result = ifelse(topology_results_selected_analytical_df[[selected_fill_parameter]] == selected_parameter, topology_results_selected_analytical_df$model_result / dec_smax_model$L, topology_results_selected_analytical_df$model_result)
    #         predicted_results_df$model_result = ifelse(predicted_results_df$fill_parameter == selected_parameter, predicted_results_df$model_result / dec_smax_model$L, predicted_results_df$model_result)
    #         results_selected_df[[input$boxplot_parameter]] = ifelse(results_selected_df[[selected_fill_parameter]] == selected_parameter, results_selected_df[[input$boxplot_parameter]] / dec_smax_model$L, results_selected_df[[input$boxplot_parameter]])
    #         topology_results_selected_by_size_df$mean_norm = ifelse(topology_results_selected_by_size_df[[selected_fill_parameter]] == selected_parameter, topology_results_selected_by_size_df$mean_norm / dec_smax_model$L, topology_results_selected_by_size_df$mean_norm)
    #       }
    #     }
    #   }
    #   predicted_results_df = predicted_results_df %>% rename(!!input$boxplot_parameter := model_result)
    #   if(("fill_parameter" %in% colnames(predicted_results_df)) && (!(is.null(selected_fill_parameter))) && (!(selected_fill_parameter == "type_dataset"))){
    #     predicted_results_df = predicted_results_df %>% rename(!!selected_fill_parameter := fill_parameter)
    #   }
    #   
    #   # Calculate relative error
    #   topology_results_selected_pred_vs_reality_df = topology_results_selected_analytical_df %>% inner_join(topology_results_selected_by_size_df %>% ungroup() %>% select("size", !!combination_metric, !!selected_fill_parameter), by=c("size", selected_fill_parameter)) %>% unique()
    #   topology_results_selected_pred_vs_reality_df$relative.error = abs((topology_results_selected_pred_vs_reality_df[[combination_metric]] - topology_results_selected_pred_vs_reality_df$model_result)) / topology_results_selected_pred_vs_reality_df[[combination_metric]]
    #   # Prepare output table
    #   analytical_model_summary_df = topology_results_selected_pred_vs_reality_df %>% 
    #     filter(model %in% c(input$topology_type_analytical_model)) %>% 
    #     select(one_of(selected_fill_parameter), "model", "formula", "adj.r.squared", "relative.error") %>% 
    #     unique() %>%
    #     # Calculate relative error mean
    #     group_by_at(c(selected_fill_parameter, "model", "formula", "adj.r.squared")) %>%
    #     mutate(relative.error.mean=mean(relative.error)) %>%
    #     select(!("relative.error")) %>%
    #     unique() %>%
    #     ungroup()
    #   # Remove table if it is there
    #   #removeUI(selector = "#topologyAnalyticalModelTableID")
    #   # Insert updated table
    #   #insertUI(
    #   #  selector = '#topologyAnalyticalModelTable',
    #   #  # wrap element in a div with id for ease of removal
    #   #  ui = tags$div(
    #   #    renderTable(analytical_model_summary_df), 
    #   #    id = "topologyAnalyticalModelTableID"
    #   #  )
    #   #)
    #   # Render table
    #   output$topologyAnalyticalModelTable = renderTable(analytical_model_summary_df)
    # } else {
    #   # Remove table if it is rendered
    #   #output$topologyAnalyticalModelTable <- renderUI({return(NULL)})
    #   removeUI(selector = "#topologyAnalyticalModelTableID")
    #   # Remove table if it is there
    #   #removeUI(selector = "#topologyAnalyticalModelTableID")
    # }
    
    # Cut values before
    #results_selected_df = results_selected_df %>% filter(size <= 2500)
    #topology_results_selected_by_size_df = topology_results_selected_by_size_df %>% filter(size <= 2500)
    #topology_results_selected_analytical_df = topology_results_selected_analytical_df %>% filter(size <= 2500)
    #predicted_results_df = predicted_results_df %>% filter(size <= 2500)
    
    return(list(results_selected_df=results_selected_df, 
                topology_results_selected_by_size_df=topology_results_selected_by_size_df, 
                topology_results_selected_analytical_df=topology_results_selected_analytical_df, 
                stretched_exponential_regression_df=stretched_exponential_regression_df,
                predicted_results_df=predicted_results_df,
                analytical_model_summary_df = analytical_model_summary_df,
                selected_fill_parameter=selected_fill_parameter, 
                combination_metric=combination_metric))
  })

  output$topologyBoxPlot <- renderPlot({
    
    processed_inputs = process_inputs()
    results_selected_df=processed_inputs$results_selected_df
    topology_results_selected_by_size_df=processed_inputs$topology_results_selected_by_size_df
    topology_results_selected_analytical_df = processed_inputs$topology_results_selected_analytical_df %>% filter(model %in% input$topology_type_analytical_model)
    stretched_exponential_regression_df = processed_inputs$stretched_exponential_regression_df  %>% filter(model %in% input$topology_type_analytical_model)
    predicted_results_df = processed_inputs$predicted_results_df  %>% filter(model %in% input$topology_type_analytical_model)
    analytical_model_summary_df = processed_inputs$analytical_model_summary_df
    selected_fill_parameter=processed_inputs$selected_fill_parameter
    combination_metric=processed_inputs$combination_metric

    #------------------#
    # Create main plot #
    #------------------#
    
    if ((!(isTruthy(input$topology_analytical))) | ((isTruthy(input$topology_analytical)) & (input$topology_type_analytical_model_output == "fit.plot"))){
      
      #----------------------------------#
      # Define main elements of the plot #
      #----------------------------------#
      
      if (is.null(selected_fill_parameter)){
        # Plot without fill
        topology_results_plot = ggplot(results_selected_df, aes(x=size, y=.data[[input$boxplot_parameter]])) + 
          guides(fill=guide_legend(title="")) +
          geom_point(alpha=0.5, size=3, col=2, fill=2) +
          geom_line(data = topology_results_selected_by_size_df,
                    aes(x = size, y = get(combination_metric), group=1),
                    col=2, lwd=1) +
          geom_line(data = topology_results_selected_analytical_df,
                    aes(x = size, y = model_result, group=model, col=model),
                    lwd=1) +
          scale_color_brewer(name = "model", palette = "Accent", direction = 1)
          #scale_color_manual(name = "model", values = c("#1f78b4", "#33a02c"))
  
      } else {
        # Plot with fill
        topology_results_plot = ggplot(results_selected_df, aes(x=size, y=.data[[input$boxplot_parameter]], col=get(selected_fill_parameter))) + 
          guides(col=guide_legend(title=parameter2label[[selected_fill_parameter]])) +
          geom_point(alpha=0.5, size=3) +
          geom_line(data = topology_results_selected_by_size_df,
                    aes(x = size, y = get(combination_metric)),
                    lwd=1) +
          geom_line(data = topology_results_selected_analytical_df %>% filter(model %in% c(input$topology_type_analytical_model)),
                    aes(x = size, y = model_result, group=get(selected_fill_parameter)),
                    lwd=1, col="black")
        
      }
      
      #----------------#
      # Customize plot #
      #----------------#
      
      # Define labels
      label_x = "Number of samples"
      label_y = parameter2label[[input$boxplot_parameter]]
      

      # Change axes label if re-scaling
      if(isTRUE(input$topology_normalize_x)){
        label_x= paste(label_x, " (Norm.)", sep="")
      }
      if(isTRUE(input$topology_normalize_y)){
        if(input$topology_type_normalization == "divide.max.value"){
          label_y= paste(label_y, " (Norm.)", sep="")
        } else if (input$topology_type_normalization == "divide.max.possible.value"){
          label_y= paste(label_y, " / Max. num. edges", sep="")
          #label_y= paste(label_y, " / (V*(V-1)/2)", sep="")
        } else if(input$topology_type_normalization == "divide.L") {
          label_y= paste(label_y, " / L", sep="")
        }
      }
      
      # Plot p-value threshold line (if a parameter involving p-value is selected)
      if (input$boxplot_parameter == "log_disease_lcc_pvalue"){
        topology_results_plot = topology_results_plot +
          geom_hline(yintercept=abs(log10(as.numeric(input$topology_pvalue_threshold))), linetype="dashed", color="red", size=1)
      }
      if (input$boxplot_parameter == "disease_lcc_pvalue"){
        topology_results_plot = topology_results_plot +
          geom_hline(yintercept=as.numeric(input$topology_pvalue_threshold), linetype="dashed", color="red", size=1)
      }
      
      # Plot standard deviation (if requested by the user)
      if(isTruthy(input$topology_sd)){
        topology_results_plot$layers[[1]] = NULL # Remove geom_point()
        label_y= paste(label_y, " (SD)", sep="") # Change y axis label
      }
      
      # Transform to log scale (if requested by the user)
      if (isTruthy(input$topology_log_x)){
        label_x= paste(label_x, " (Log)", sep="")
        topology_results_plot = topology_results_plot + 
          scale_x_continuous(trans = scales::log10_trans())
      }
      if (isTruthy(input$topology_log_y)){
        label_y= paste(label_y, " (Log)", sep="")
        topology_results_plot = topology_results_plot + 
          scale_y_continuous(trans = scales::log10_trans())
      }
      
      # Customize axes/labels/grid
      topology_results_plot = topology_results_plot + 
        theme_linedraw() +
        xlab(label_x) +
        ylab(label_y) +
        theme(plot.title =  element_text(size = 17, face="bold"), axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 15), legend.text = element_text(size = 14), legend.title=element_text(size=15, face="bold"), legend.position="bottom")
      
      topology_results_plot
    
    } else {
      
      if ((isTruthy(input$topology_analytical)) & (input$topology_type_analytical_model_output == "prediction.plot")){

        #------------------------#
        # Create prediction plot #
        #------------------------#
        
        # Change axes label if re-scaling
        label_x = "Number of samples"
        label_y = parameter2label[[input$boxplot_parameter]]
        
        # Change axes label if re-scaling
        if(isTRUE(input$topology_normalize_x)){
          label_x= paste(label_x, " (Norm.)", sep="")
        }
        if(isTRUE(input$topology_normalize_y)){
          if(input$topology_type_normalization == "divide.max.value"){
            label_y= paste(label_y, " (Norm.)", sep="")
          } else if (input$topology_type_normalization == "divide.max.possible.value"){
            label_y= paste(label_y, " / Max. num. edges", sep="")
            #label_y= paste(label_y, " / (V*(V-1)/2)", sep="")
          } else if(input$topology_type_normalization == "divide.L") {
            label_y= paste(label_y, " / L", sep="")
          }
        }
        
        # Plot log(gradient) vs log(sample size)
        if (is.null(selected_fill_parameter)){
          # Plot without fill
          prediction_plot = ggplot(predicted_results_df, aes(x=size, y=model_result)) +
            geom_line(size=2, col=2)
        } else {
          # Plot with fill
          prediction_plot = ggplot(predicted_results_df, aes(x=size, y=model_result, col=get(selected_fill_parameter))) +
            geom_line(size=2) +
            guides(col=guide_legend(title=parameter2label[[selected_fill_parameter]]))
        }
        
        # Transform to log scale (if requested by the user)
        if (isTruthy(input$topology_log_x)){
          label_x= paste(label_x, " (Log)", sep="")
          prediction_plot = prediction_plot + 
            scale_x_continuous(trans = scales::log10_trans())
        }
        if (isTruthy(input$topology_log_y)){
          label_y= paste(label_y, " (Log)", sep="")
          prediction_plot = prediction_plot + 
            scale_y_continuous(trans = scales::log10_trans())
        }
        
        prediction_plot = prediction_plot +
          theme_linedraw() +
          xlab(label_x) +
          ylab(label_y) +
          theme(plot.title =  element_text(size = 17, face="bold"), axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 15), legend.text = element_text(size = 14), legend.title=element_text(size=15, face="bold"), legend.position="bottom")
        
        prediction_plot

      } else if ((isTruthy(input$topology_analytical)) & ((input$topology_type_analytical_model_output == "regression.plot") | (input$topology_type_analytical_model_output == "cumulative.regression.plot"))){
        
        if(((isTruthy(input$topology_analytical)) & (!("Logarithmic" %in% c(input$topology_type_analytical_model))) & (!(length(input$topology_type_analytical_model) == 0))) | ((isTruthy(input$topology_analytical)) & ("Logarithmic" %in% c(input$topology_type_analytical_model)) & (length(c(input$topology_type_analytical_model)) > 1))){
          #if(!((isTruthy(input$topology_analytical)) & ("Logarithmic" %in% c(input$topology_type_analytical_model)) & (length(c(input$topology_type_analytical_model)) == 1))){
          
          #------------------------#
          # Create regression plot #
          #------------------------#
          
          # Plot log(gradient) vs log(sample size)
          if (is.null(selected_fill_parameter)){
            # Plot without fill
            # Make it cumulative
            if(input$topology_type_analytical_model_output == "cumulative.regression.plot"){
              regression_plot  = ggplot(stretched_exponential_regression_df, aes(x=x, y=y))
              regression_plot = regression_plot + geom_line(aes(y = 1 - ..y..), stat='ecdf', size=2, col=2) # reverse plot
            } else {
              regression_plot = ggplot(stretched_exponential_regression_df, aes(x=x, y=y)) +
                geom_point(alpha=0.9, size=3, col=2) +
                geom_line(aes(x = x, y = regression_line), col=2)
            }
          } else {
            # Plot with fill
            # Make it cumulative
            if(input$topology_type_analytical_model_output == "cumulative.regression.plot"){
              regression_plot  = ggplot(stretched_exponential_regression_df, aes(x=x, y=y, col=get(selected_fill_parameter)))
              regression_plot = regression_plot + geom_line(aes(y = 1 - ..y..), stat='ecdf', size=2) + # reverse plot
                guides(col=guide_legend(title=parameter2label[[selected_fill_parameter]]))
            } else {
              regression_plot = ggplot(stretched_exponential_regression_df, aes(x=x, y=y, col=get(selected_fill_parameter))) + 
                guides(col=guide_legend(title=parameter2label[[selected_fill_parameter]])) +
                geom_point(alpha=0.9, size=3) +
                geom_line(aes(x = x, y = regression_line, col=get(selected_fill_parameter)))
            }
          }
          
          if((input$topology_type_analytical_model == "Stretched exponential (by optimization)") | (input$topology_type_analytical_model == "Stretched exponential (by linear fit)")){
            label_y = "Ln( Ln(smax) - Ln(s) )"
          } else if(input$topology_type_analytical_model == "Stretched exponential (without L)"){
            label_y = "Ln(P(N))"
          }
          
          regression_plot = regression_plot +
            theme_linedraw() +
            xlab("Ln(N)") +
            ylab(label_y) +
            ggtitle("Linear regression of the analytical model") +
            theme(plot.title =  element_text(size = 17, face="bold"), axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 15), legend.text = element_text(size = 14), legend.title=element_text(size=15, face="bold"), legend.position="bottom")
          
          regression_plot
          
        }
        
      }
      
    }
    
  }, height = 500, width = 800)

  
  
  
  
  
  output$regressionPlot <- renderPlot({
    
    processed_inputs = process_inputs()
    stretched_exponential_regression_df=processed_inputs$stretched_exponential_regression_df
    selected_fill_parameter=processed_inputs$selected_fill_parameter
    if(((isTruthy(input$topology_analytical)) & (!("Logarithmic" %in% c(input$topology_type_analytical_model))) & (!(length(input$topology_type_analytical_model) == 0))) | ((isTruthy(input$topology_analytical)) & ("Logarithmic" %in% c(input$topology_type_analytical_model)) & (length(c(input$topology_type_analytical_model)) > 1))){
      #if(!((isTruthy(input$topology_analytical)) & ("Logarithmic" %in% c(input$topology_type_analytical_model)) & (length(c(input$topology_type_analytical_model)) == 1))){
      # Plot log(gradient) vs log(sample size)
      if (is.null(selected_fill_parameter)){
        # Plot without fill
        # Make it cumulative
        if(isTruthy(input$topology_cumulative)){
          regression_plot  = ggplot(stretched_exponential_regression_df, aes(x=x, y = y))
          regression_plot = regression_plot + geom_line(aes(y = 1 - ..y..), stat='ecdf', size=2, col=2) # reverse plot
        } else {
          regression_plot = ggplot(stretched_exponential_regression_df, aes(x=x, y=y)) +
            geom_point(alpha=0.9, size=3, col=2) +
            geom_line(aes(x = x, y = regression_line), col=2)
        }
      } else {
        # Plot with fill
        # Make it cumulative
        if(isTruthy(input$topology_cumulative)){
          regression_plot  = ggplot(stretched_exponential_regression_df, aes(x=x, y = y, col=get(selected_fill_parameter)))
          regression_plot = regression_plot + geom_line(aes(y = 1 - ..y..), stat='ecdf', size=2) + # reverse plot
            guides(col=guide_legend(title=parameter2label[[selected_fill_parameter]]))
        } else {
          regression_plot = ggplot(stretched_exponential_regression_df, aes(x=x, y=y, col=get(selected_fill_parameter))) + 
            guides(col=guide_legend(title=parameter2label[[selected_fill_parameter]])) +
            geom_point(alpha=0.9, size=3) +
            geom_line(aes(x = x, y = regression_line, col=get(selected_fill_parameter)))
        }
      }
      
      if((input$topology_type_analytical_model == "Stretched exponential (by optimization)") | (input$topology_type_analytical_model == "Stretched exponential (by linear fit)")){
        label_y = "Ln( Ln(smax) - Ln(s) )"
      } else if(input$topology_type_analytical_model == "Stretched exponential (without L)"){
        label_y = "Ln(P(N))"
      }
      
      regression_plot = regression_plot +
        theme_linedraw() +
        xlab("Ln(N)") +
        ylab(label_y) +
        ggtitle("Linear regression of the analytical model") +
        theme(plot.title =  element_text(size = 17, face="bold"), axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 15), legend.text = element_text(size = 14), legend.title=element_text(size=15, face="bold"), legend.position="bottom")
      
      regression_plot
      
    }
  })
  
  # conditional UI: we only create the space for the plot when the conditions are met
  output$regressionPlotID = renderUI({
    if(((isTruthy(input$topology_analytical)) & (!("Logarithmic" %in% c(input$topology_type_analytical_model))) & (!(length(input$topology_type_analytical_model) == 0))) | ((isTruthy(input$topology_analytical)) & ("Logarithmic" %in% c(input$topology_type_analytical_model)) & (length(c(input$topology_type_analytical_model)) > 1))){
      plotOutput(outputId = "regressionPlot", height = 500, width = 800)
    }
  })
  
  output$predPlot <- renderPlot({
    
    processed_inputs = process_inputs()
    predicted_results_df=processed_inputs$predicted_results_df %>% filter(model %in% c(input$topology_type_analytical_model))
    selected_fill_parameter=processed_inputs$selected_fill_parameter
    
    # Change axes label if re-scaling
    label_x = "Number of samples"
    label_y = parameter2label[[input$boxplot_parameter]]
    if(isTRUE(input$topology_normalize_x)){
      label_x= paste(label_x, " (Norm.)", sep="")
    }
    if(isTRUE(input$topology_normalize_y)){
      label_y= paste(label_y, " (Norm.)", sep="")
    }
    
    if(((isTruthy(input$topology_analytical)) & (!("Logarithmic" %in% c(input$topology_type_analytical_model))) & (!(length(input$topology_type_analytical_model) == 0))) | ((isTruthy(input$topology_analytical)) & ("Logarithmic" %in% c(input$topology_type_analytical_model)) & (length(c(input$topology_type_analytical_model)) > 1))){
      # Plot log(gradient) vs log(sample size)
      if (is.null(selected_fill_parameter)){
        # Plot without fill
        prediction_plot = ggplot(predicted_results_df, aes(x=size, y=.data[[input$boxplot_parameter]])) +
          geom_line(size=2, col=2)
      } else {
        # Plot with fill
        prediction_plot = ggplot(predicted_results_df, aes(x=size, y=.data[[input$boxplot_parameter]], col=get(selected_fill_parameter))) +
          geom_line(size=2) +
          guides(col=guide_legend(title=parameter2label[[selected_fill_parameter]]))
      }
      
      # Transform to log scale (if requested by the user)
      if (isTruthy(input$topology_log_x)){
        label_x= paste(label_x, " (Log)", sep="")
        prediction_plot = prediction_plot + 
          scale_x_continuous(trans = scales::log10_trans())
      }
      if (isTruthy(input$topology_log_y)){
        label_y= paste(label_y, " (Log)", sep="")
        prediction_plot = prediction_plot + 
          scale_y_continuous(trans = scales::log10_trans())
      }
      
      prediction_plot = prediction_plot +
        theme_linedraw() +
        xlab(label_x) +
        ylab(label_y) +
        theme(plot.title =  element_text(size = 17, face="bold"), axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 15), legend.text = element_text(size = 14), legend.title=element_text(size=15, face="bold"), legend.position="bottom")
      
      prediction_plot
      
    }
  })
  
  # conditional UI: we only create the space for the plot when the conditions are met
  output$predPlotID = renderUI({
    if(((isTruthy(input$topology_analytical)) & (!("Logarithmic" %in% c(input$topology_type_analytical_model))) & (!(length(input$topology_type_analytical_model) == 0))) | ((isTruthy(input$topology_analytical)) & ("Logarithmic" %in% c(input$topology_type_analytical_model)) & (length(c(input$topology_type_analytical_model)) > 1))){
      plotOutput(outputId = "predPlot", height = 500, width = 800)
    }
  })
  
}

