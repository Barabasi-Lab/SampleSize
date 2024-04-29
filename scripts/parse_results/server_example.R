library(data.table)
library(dplyr)
library(ggplot2)
library(stats)

# Input files
setwd("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeShiny")
input_dir = 'data'
topology_results_file = paste(input_dir, 'analysis_topology.csv', sep='/')
ppi_results_file = paste(input_dir, 'analysis_ppi.csv', sep='/')
disease_genes_results_file = paste(input_dir, 'analysis_disease_genes.csv', sep='/')
essential_genes_results_file = paste(input_dir, 'analysis_essential_genes.csv', sep='/')
numbers_complete_graph_file = paste(input_dir, 'dataset_numbers_complete_graph.txt', sep='/')
methods_selected <- c("pearson", "aracne")
thresholds_selected <- c(0.05, 0)
type_data_selection = 'pearson_aracne'
#methods_selected <- c("pearson")
#thresholds_selected <- c(0.05)
#type_data_selection = 'pearson_pval_0.05'
output_dir = paste('/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeShiny/data/example_', type_data_selection, sep = '')
dir.create(output_dir, showWarnings = FALSE)

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
  # Calculate mean of repetitions from same sample size
  results_dataframe = results_dataframe %>% 
    arrange(get(x_parameter), rep) %>%
    group_by(get(x_parameter)) %>%
    summarise_at(vars(all_of(y_parameter)), list(mean=mean, median=median, sd=sd)) %>%
    rename(!!x_parameter := "get(x_parameter)") %>%
    rename(!!y_parameter := "mean") %>%
    ungroup()
  # Calculate linear fit
  lm_summary = summary(lm(get(y_parameter)~log(get(x_parameter)), data=results_dataframe))
  used_data=data.frame(y=results_dataframe[[y_parameter]], x=log(results_dataframe[[x_parameter]]))
  return(list(lm_summary=lm_summary, used_data=used_data, slope=coef(lm_summary)[2], intercept=coef(lm_summary)[1], a=coef(lm_summary)[2], b=coef(lm_summary)[1], L=NaN, adj.r.squared=lm_summary$adj.r.squared))
}

#'  get_exponential_tendency
#'  Method to obtain a logarithmic fit from the data.
#'  @param results_dataframe Dataframe containing the data.
#'  @param y_parameter Parameter of the Y axis.
#'  @param x_parameter Parameter of the X axis.
#'  
get_exponential_tendency = function(results_dataframe, y_parameter, x_parameter){
  # Calculate mean of repetitions from same sample size
  results_dataframe = results_dataframe %>% 
    arrange(get(x_parameter), rep) %>%
    group_by(get(x_parameter)) %>%
    summarise_at(vars(all_of(y_parameter)), list(mean=mean, median=median, sd=sd)) %>%
    rename(!!x_parameter := "get(x_parameter)") %>%
    rename(!!y_parameter := "mean") %>%
    ungroup()
  # Calculate linear fit
  lm_summary = summary(lm(log(get(y_parameter))~get(x_parameter), data=results_dataframe))
  used_data=data.frame(y=log(results_dataframe[[y_parameter]]), x=results_dataframe[[x_parameter]])
  return(list(lm_summary=lm_summary, used_data=used_data, slope=coef(lm_summary)[2], intercept=coef(lm_summary)[1], a=coef(lm_summary)[2], b=coef(lm_summary)[1], L=NaN, adj.r.squared=lm_summary$adj.r.squared))
}

#'  get_linear_tendency
#'  Method to obtain a logarithmic fit from the data.
#'  @param results_dataframe Dataframe containing the data.
#'  @param y_parameter Parameter of the Y axis.
#'  @param x_parameter Parameter of the X axis.
#'  
get_linear_tendency = function(results_dataframe, y_parameter, x_parameter){
  # Calculate mean of repetitions from same sample size
  results_dataframe = results_dataframe %>% 
    arrange(get(x_parameter), rep) %>%
    group_by(get(x_parameter)) %>%
    summarise_at(vars(all_of(y_parameter)), list(mean=mean, median=median, sd=sd)) %>%
    rename(!!x_parameter := "get(x_parameter)") %>%
    rename(!!y_parameter := "mean") %>%
    ungroup()
  # Calculate linear fit
  lm_summary = summary(lm(get(y_parameter)~get(x_parameter), data=results_dataframe))
  used_data=data.frame(y=results_dataframe[[y_parameter]], x=results_dataframe[[x_parameter]])
  return(list(lm_summary=lm_summary, used_data=used_data, slope=coef(lm_summary)[2], intercept=coef(lm_summary)[1], a=coef(lm_summary)[2], b=coef(lm_summary)[1], L=NaN, adj.r.squared=lm_summary$adj.r.squared))
}

#'  get_square_root_tendency
#'  Method to obtain a logarithmic fit from the data.
#'  @param results_dataframe Dataframe containing the data.
#'  @param y_parameter Parameter of the Y axis.
#'  @param x_parameter Parameter of the X axis.
#'  
get_square_root_tendency = function(results_dataframe, y_parameter, x_parameter){
  # Calculate mean of repetitions from same sample size
  results_dataframe = results_dataframe %>% 
    arrange(get(x_parameter), rep) %>%
    group_by(get(x_parameter)) %>%
    summarise_at(vars(all_of(y_parameter)), list(mean=mean, median=median, sd=sd)) %>%
    rename(!!x_parameter := "get(x_parameter)") %>%
    rename(!!y_parameter := "mean") %>%
    ungroup()
  # Calculate linear fit
  lm_summary = summary(lm(log(get(y_parameter))~log(get(x_parameter)), data=results_dataframe))
  used_data=data.frame(y=log(results_dataframe[[y_parameter]]), x=log(results_dataframe[[x_parameter]]))
  return(list(lm_summary=lm_summary, used_data=used_data, slope=coef(lm_summary)[2], intercept=coef(lm_summary)[1], a=coef(lm_summary)[2], b=coef(lm_summary)[1], L=NaN, adj.r.squared=lm_summary$adj.r.squared))
}

#'  get_exponential_decay_tendency
#'  Method to obtain an exponential decay fit from the data.
#'  @param results_dataframe Dataframe containing the data.
#'  @param y_parameter Parameter of the Y axis.
#'  @param x_parameter Parameter of the X axis.
#'  
get_exponential_decay_tendency = function(results_dataframe, y_parameter, x_parameter, L){
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
  # Calculate polynomial equation = gradient vs. sample size (starting from 2nd position) (removing negatives)
  lm_summary = summary(lm(frac_grad[boo]~results_dataframe_mean[[x_parameter]][2:length(results_dataframe_mean[[x_parameter]])][boo]))
  used_data = data.frame(y=frac_grad[boo], x=results_dataframe_mean[[x_parameter]][2:length(results_dataframe_mean[[x_parameter]])][boo])
  return(list(lm_summary=lm_summary, used_data=used_data, slope=coef(lm_summary)[2], intercept=coef(lm_summary)[1], a=coef(lm_summary)[2], b=coef(lm_summary)[1], L=L, adj.r.squared=lm_summary$adj.r.squared))
}

#'  get_gaussian_tendency
#'  Method to obtain a gaussian fit from the data.
#'  @param results_dataframe Dataframe containing the data.
#'  @param y_parameter Parameter of the Y axis.
#'  @param x_parameter Parameter of the X axis.
#'  
get_gaussian_tendency = function(results_dataframe, y_parameter, x_parameter){
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
  # Calculate polynomial equation = gradient vs. sample size (starting from 2nd position) (removing negatives)
  lm_summary = summary(lm(frac_grad[boo]~results_dataframe_mean[[x_parameter]][2:length(results_dataframe_mean[[x_parameter]])][boo]))
  used_data = data.frame(y=frac_grad[boo], x=results_dataframe_mean[[x_parameter]][2:length(results_dataframe_mean[[x_parameter]])][boo])
  return(list(lm_summary=lm_summary, used_data=used_data, slope=coef(lm_summary)[2], intercept=coef(lm_summary)[1], a=coef(lm_summary)[2], b=coef(lm_summary)[1], L=NaN, adj.r.squared=lm_summary$adj.r.squared))
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
  return(list(lm_summary=lm_summary, used_data=used_data, slope=coef(lm_summary)[2], intercept=coef(lm_summary)[1], a=coef(lm_summary)[2], b=coef(lm_summary)[1], L=NaN, adj.r.squared=lm_summary$adj.r.squared))
}

#'  calculate_predictions_using_stretched_exponential_model_without_L
#'  Formula to calculate stretched exponential of significant interactions from a list of sample sizes without using L.
#'  This formula does not require the use of L
#'  @param x List of sample sizes.
#'  @param a
#'  @param b
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
  return(list(lm_summary=lm_summary, used_data=used_data, slope=coef(lm_summary)[2], intercept=coef(lm_summary)[1], a=a, b=b, L=L, adj.r.squared=lm_summary$adj.r.squared))
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
  #return(list(lm_summary=lm_summary, used_data=used_data, a=coef(lm_summary)[2], b=coef(lm_summary)[1], L=L, adj.r.squared=lm_summary$adj.r.squared))
  return(list(lm_summary=lm_summary, used_data=used_data, slope=coef(lm_summary)[2], intercept=coef(lm_summary)[1], a=a, b=b, L=L, adj.r.squared=lm_summary$adj.r.squared))
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

#'  calculate_predictions_using_stretched_exponential_model_optimized
#'  Formula to calculate exponential decay of significant interactions from a list of sample sizes.
#'  This formula requires the use of the L parameter
#'  @param x List of sample sizes.
#'  @param a Slope coefficient.
#'  @param b Intercept coefficient.
#'  
calculate_predictions_using_exponential_decay_model = function(x, L, a, b){
  y = L * exp(-a * x)
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
  } else if(model == "Linear"){
    model_output = get_linear_tendency(results_dataframe=results_dataframe, y_parameter=y_parameter, x_parameter=x_parameter)
    model_result = (sort(unique(results_dataframe[[x_parameter]]))*model_output$a + model_output$b)
  } else if(model == "Exponential"){
    model_output = get_exponential_tendency(results_dataframe=results_dataframe, y_parameter=y_parameter, x_parameter=x_parameter)
    model_result = (exp(sort(unique(results_dataframe[[x_parameter]]))*model_output$a + model_output$b))
  } else if(model == "Square root"){
    model_output = get_square_root_tendency(results_dataframe=results_dataframe, y_parameter=y_parameter, x_parameter=x_parameter)
    model_result = (exp((log(sort(unique(results_dataframe[[x_parameter]])))*model_output$a + model_output$b)))
  } else if(model == "Exponential decay"){
    model_output = get_exponential_decay_tendency(results_dataframe=results_dataframe, y_parameter=y_parameter, x_parameter=x_parameter, L=L)
    model_result = calculate_predictions_using_exponential_decay_model(x=sort(unique(results_dataframe[[x_parameter]])), L=L, a=model_output$a, b=model_output$b)
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
  } else if(model == "Linear"){
    prediction_result = (x_list*a + b)
  } else if(model == "Exponential"){
    prediction_result = (exp((x_list*a + b)))
  } else if(model == "Square root"){
    prediction_result = (exp((log(x_list)*a + b)))
  } else if((model == "Stretched exponential (by optimization)") | (model == "Stretched exponential (by linear fit)")){
    prediction_result = calculate_predictions_using_stretched_exponential_model_optimized(x=x_list, L=L, a=a, b=b)
  } else if(model == "Stretched exponential (without L)"){
    prediction_result = calculate_predictions_using_stretched_exponential_model_without_L(x=x_list, a=a, b=b)
  } else if(model == "Exponential decay"){
    prediction_result = calculate_predictions_using_exponential_decay_model(x=x_list, L=L, a=a, b=b)
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
  } else if(model == "Linear"){
    formula_name = paste("F(x) = ", formatC(round(a, 2), format = "e", digits = 2), " * x + ", formatC(round(b, 2), format = "e", digits = 2), sep="")
  } else if(model == "Exponential"){
    formula_name = paste("F(x) = exp(", formatC(round(a, 2), format = "e", digits = 2), " * x + ", formatC(round(b, 2), format = "e", digits = 2), ")", sep="")
  } else if(model == "Square root"){
    formula_name = paste("F(x) = exp(", formatC(round(a, 2), format = "e", digits = 2), " * ln(x) + ", formatC(round(b, 2), format = "e", digits = 2), ")", sep="")
  } else if((model == "Stretched exponential (by optimization)") | (model == "Stretched exponential (by linear fit)")){
    formula_name = paste("F(x) = ", round(L, 2), " * exp[(", round(b, 2), " * x**((", round(a, 2), ") + 1)) / ((", round(a, 2), ") + 1) ]", sep="")
  } else if(model == "Stretched exponential (without L)"){
    formula_name = paste("F(x) = exp[ (exp(", round(b, 2), ") * x**(1 + (", round(a, 2), "))) / (1 + (", round(a, 2), ")) ]", sep="")
  } else if(model == "Exponential decay"){
    formula_name = paste("F(x) = exp(-", round(a, 2), " * x)", sep="")
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
results_df = inner_join(topology_results_df, ppi_results_df, by = c("method", "dataset", "type_dataset", "subclassification", "size", "rep", "type_correlation", "threshold")) %>% inner_join(disease_genes_results_df, by = c("method", "dataset", "type_dataset", "subclassification", "size", "rep", "type_correlation", "threshold")) %>% inner_join(essential_genes_results_df, by = c("method", "dataset", "type_dataset", "subclassification", "size", "rep", "type_correlation", "threshold"))
results_df$type_dataset = paste(results_df$dataset, results_df$type_dataset, sep=":") # Join dataset and type_dataset
results_df$type_dataset = tolower(results_df$type_dataset)
results_df$type_dataset = ifelse(results_df$subclassification == "normal", paste(results_df$type_dataset, "normal", sep="-"), results_df$type_dataset)
results_df$threshold = as.character(results_df$threshold) # Consider threshold as a discrete variable
results_df$log_disease_lcc_pvalue = abs(log10(results_df$disease_lcc_pvalue))
results_df$log_disease_lcc_pvalue = replace(results_df$log_disease_lcc_pvalue, is.infinite(results_df$log_disease_lcc_pvalue),abs(log(0.000001))) # Replace infinite values by very high values
numbers_complete_graph_df = fread(numbers_complete_graph_file) %>% rename("type_dataset" = "dataset")
numbers_complete_graph_df$type_dataset = tolower(numbers_complete_graph_df$type_dataset)


#--------------------#
# Select information #
#--------------------#

#results_selected_df = results_df %>% filter((type_dataset %in% c("tcga:tcga-brca")) & (method %in% c("pearson")) & (type_correlation %in% c("all")) & (threshold %in% c(0.05))) # For a test
#results_selected_df = results_df %>% filter((type_dataset %in% c("tcga:tcga-brca", "tcga:tcga-ucec")) & (method %in% c("pearson")) & (type_correlation %in% c("all")) & (threshold %in% c(0.05))) # For a test with multiple parameters
#results_selected_df = results_df %>% filter((method %in% c("pearson")) & (type_correlation %in% c("all")) & (threshold %in% c(0.05))) # For a test with multiple parameters
results_selected_df = results_df %>% filter((method %in% methods_selected) & (type_correlation %in% c("all")) & (threshold %in% thresholds_selected)) # Include the methods and thresholds selected by the user

# Merge with total number of edges
results_selected_df = results_selected_df %>% inner_join(numbers_complete_graph_df %>% select("type_dataset", "total_num_edges"), by = "type_dataset") %>% group_by_at(c("dataset", "type_dataset")) %>% mutate(total_num_edges_norm=total_num_edges/max(num_edges)) %>% ungroup() %>% as.data.frame()

# Select by diseases
group_vars = c("type_dataset", "method", "size", "type_correlation", "threshold", "total_num_edges")
type_analysis = "topology"
if(type_analysis == "disease_genes"){
  diseases_selected <- c("alzheimer disease", "arthritis rheumatoid", "cardiomyopathies", "diabetes mellitus type 2")
  results_selected_df %<>% filter(disease %in% diseases_selected) %>% select(-disease_class) %>% unique()
  group_vars = c(group_vars, "disease")
} else {
  results_selected_df %<>% select(!(c("num_disease_genes", "num_disease_edges", "fraction_disease_genes", "num_disease_components", "disease", "disease_class", "num_disease_lcc_nodes", "num_disease_lcc_edges", "fraction_disease_lcc_nodes", "disease_lcc_z", "disease_lcc_pvalue", "log_disease_lcc_pvalue"))) %>% unique()
}


#--------------------------#
# Calculate mean / sd line #
#--------------------------#

# Calculate mean and standard deviation tendency
boxplot_parameter = "num_edges"
topology_results_selected_by_size_df = results_selected_df %>%
  group_by_at(group_vars) %>%
  summarise_at(vars(all_of(boxplot_parameter)), list(mean=mean, sd=sd)) %>%
  arrange(size)

# Check if the user wants to plot mean or standard deviation
combination_metric = "mean"

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
cols = c("model", "model_result", "size", "slope", "intercept", "a", "b", "L", "adj.r.squared", "max_value_in_dataset", "formula", "type_dataset", "fill_parameter")
topology_results_selected_analytical_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
cols = c("model", "y", "x", "regression", "type_dataset", "fill_parameter")
stretched_exponential_regression_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
cols = c("model", "model_result", "size", "max_value_in_dataset", "type_dataset", "fill_parameter")
predicted_results_df = data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols)))
analytical_model_summary_df = data.frame()
N_vals = seq(10, 50000, 10)

types_analytical_model = c("Stretched exponential (by optimization)", "Stretched exponential (by linear fit)", "Stretched exponential (without L)", "Logarithmic", "Exponential", "Linear", "Square root", "Exponential decay")
for(model in types_analytical_model){
  if (is.null(selected_fill_parameter)){
    # Calculate the analytical model
    output_variable = calculate_analytical_model(results_dataframe=results_selected_df, y_parameter=boxplot_parameter, x_parameter="size", model=model, L=unique(results_selected_df$total_num_edges_norm))
    # Calculate prediction
    prediction_result = calculate_prediction_from_analytical_model(model=model, x_list=N_vals, a=output_variable$model_output$a, b=output_variable$model_output$b, L=output_variable$model_output$L)
    # Un-normalize y axis of analytical curve (if necessary), because the analytical curve of stretched exponential is normalized by default
    max_value_in_dataset = max(results_selected_df[[boxplot_parameter]])
    if (!(model == "Logarithmic")){
      output_variable$model_result = output_variable$model_result * max_value_in_dataset
      prediction_result = prediction_result * max_value_in_dataset
      # Save data used to create the model
      stretched_exponential_regression_df = rbind(stretched_exponential_regression_df, cbind(data.frame(model=model), output_variable$model_output$used_data, data.frame(regression_line=(output_variable$model_output$slope * output_variable$model_output$used_data$x + output_variable$model_output$intercept), type_dataset=unique(results_selected_df$type_dataset)[1])))
    }
    # Calculate formula
    formula = get_formula(model=model, a=output_variable$model_output$a, b=output_variable$model_output$b, L=output_variable$model_output$L)
    # Add results
    topology_results_selected_analytical_df = rbind(topology_results_selected_analytical_df, data.frame(model=model, model_result=output_variable$model_result, size=sort(unique(results_selected_df$size)), slope=output_variable$model_output$slope, intercept=output_variable$model_output$intercept, a=output_variable$model_output$a, b=output_variable$model_output$b, L=output_variable$model_output$L, adj.r.squared=output_variable$model_output$lm_summary$adj.r.squared, max_value_in_dataset=max_value_in_dataset, formula=formula, type_dataset=unique(results_selected_df$type_dataset)[1]))
    predicted_results_df = rbind(predicted_results_df, data.frame(model=model, model_result=prediction_result, size=N_vals, max_value_in_dataset=max_value_in_dataset, type_dataset=unique(results_selected_df$type_dataset)[1]))
  } else{
    for(selected_parameter in unique(results_selected_df[[selected_fill_parameter]])){
      # Select results for a specific parameter
      results_selected_by_parameter_df = results_selected_df %>% filter(!!as.symbol(selected_fill_parameter) == selected_parameter)
      # Calculate the analytical model
      output_variable = calculate_analytical_model(results_dataframe=results_selected_by_parameter_df, y_parameter=boxplot_parameter, x_parameter="size", model=model, L=unique(results_selected_by_parameter_df$total_num_edges_norm))
      # Calculate prediction
      prediction_result = calculate_prediction_from_analytical_model(model=model, x_list=N_vals, a=output_variable$model_output$a, b=output_variable$model_output$b, L=output_variable$model_output$L)
      # Un-normalize y axis of analytical curve (if necessary), because the analytical curve of stretched exponential is normalized by default
      max_value_in_dataset = max(results_selected_by_parameter_df[[boxplot_parameter]])
      if (!(model == "Logarithmic")){
        output_variable$model_result = output_variable$model_result * max_value_in_dataset
        prediction_result = prediction_result * max_value_in_dataset
        # Save data used to create the model
        stretched_exponential_regression_df = rbind(stretched_exponential_regression_df, cbind(data.frame(model=model), output_variable$model_output$used_data, data.frame(regression_line=(output_variable$model_output$slope * output_variable$model_output$used_data$x + output_variable$model_output$intercept), type_dataset=unique(results_selected_by_parameter_df$type_dataset)[1], fill_parameter=selected_parameter)))
      }
      # Calculate formula
      formula = get_formula(model=model, a=output_variable$model_output$a, b=output_variable$model_output$b, L=output_variable$model_output$L)
      # Add results
      topology_results_selected_analytical_df = rbind(topology_results_selected_analytical_df, data.frame(model=model, model_result=output_variable$model_result, size=sort(unique(results_selected_by_parameter_df$size)), slope=output_variable$model_output$slope, intercept=output_variable$model_output$intercept, a=output_variable$model_output$a, b=output_variable$model_output$b, L=output_variable$model_output$L, adj.r.squared=output_variable$model_output$lm_summary$adj.r.squared, max_value_in_dataset=max_value_in_dataset, formula=formula, type_dataset=unique(results_selected_by_parameter_df$type_dataset), fill_parameter=selected_parameter))
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
  #filter(model %in% c(topology_type_analytical_model)) %>% 
  select(one_of(selected_fill_parameter), "type_dataset", "model", "a", "b", "L", "max_value_in_dataset", "formula", "adj.r.squared", "relative.error") %>% 
  unique() %>%
  # Calculate relative error mean
  group_by_at(c(selected_fill_parameter, "model", "formula", "adj.r.squared")) %>%
  mutate(relative.error.mean=mean(relative.error)) %>%
  select(!("relative.error")) %>%
  unique() %>%
  ungroup()
analytical_model_summary_df$unnorm_L = analytical_model_summary_df$L * analytical_model_summary_df$max_value_in_dataset
analytical_model_summary_df = analytical_model_summary_df %>% inner_join(numbers_complete_graph_df %>% select("type_dataset", "total_num_edges"), by = "type_dataset") %>% group_by_at(selected_fill_parameter) %>% mutate(density = unnorm_L/total_num_edges)


#---------------#
# Write results #
#---------------#

topology_results_file = paste(output_dir, 'topology_results_pearson_pval_0.05.txt', sep='/')
results_selected_df %>% fwrite(topology_results_file)
topology_results_by_size_file = paste(output_dir, 'topology_results_mean_pearson_pval_0.05.txt', sep='/')
topology_results_selected_by_size_df %>% fwrite(topology_results_by_size_file)
predictions_file = paste(output_dir, 'predictions_pearson_pval_0.05.txt', sep='/')
predicted_results_df %>% fwrite(predictions_file)
analytical_results_file = paste(output_dir, 'analytical_model_results_pearson_pval_0.05.txt', sep='/')
topology_results_selected_analytical_df %>% fwrite(analytical_results_file)
analytical_summary_file = paste(output_dir, 'analytical_model_summary_pearson_pval_0.05.txt', sep='/')
analytical_model_summary_df %>% fwrite(analytical_summary_file)
analytical_regression_results_file = paste(output_dir, 'analytical_model_regression_results_pearson_pval_0.05.txt', sep='/')
stretched_exponential_regression_df %>% fwrite(analytical_regression_results_file)


#------------------------------------#
# Normalize data if required by user #
#------------------------------------#

# Re-scale X axis
topology_normalize_x = FALSE
if (isTRUE(topology_normalize_x)){
  results_selected_df = results_selected_df %>% group_by_at(selected_fill_parameter) %>% mutate(max_size = max(size)) %>% mutate(norm = size/max_size) %>% select(!(all_of(c("size", "max_size")))) %>% rename(size = norm)
  topology_results_selected_by_size_df = topology_results_selected_by_size_df %>% group_by_at(selected_fill_parameter) %>% mutate(max_size=max(size)) %>% mutate(norm = size/max_size) %>% select(!(all_of(c("size", "max_size")))) %>% rename(size = norm)
  topology_results_selected_analytical_df = topology_results_selected_analytical_df %>% group_by_at(selected_fill_parameter) %>% mutate(max_size = max(size)) %>% mutate(norm = size/max_size) %>% select(!(all_of(c("size", "max_size")))) %>% rename(size = norm)
  predicted_results_df = predicted_results_df %>% group_by_at(selected_fill_parameter) %>% mutate(max_size = max(size)) %>% mutate(norm = size/max_size) %>% select(!(all_of(c("size", "max_size")))) %>% rename(size = norm)
}

# Rescale Y axis
topology_normalize_y = TRUE
topology_type_normalization = "divide.L"
topology_type_analytical_model = "Stretched exponential (by linear fit)"
if(isTRUE(topology_normalize_y)){
  if(topology_type_normalization == "divide.max.value"){
    results_selected_norm_df = results_selected_df %>% group_by_at(selected_fill_parameter) %>% mutate(max_parameter = max(get(boxplot_parameter))) %>% mutate(norm = get(boxplot_parameter)/max_parameter) %>% rename(unnorm = boxplot_parameter) %>% rename(!!boxplot_parameter := norm)
    topology_results_selected_by_size_norm_df = topology_results_selected_by_size_df %>% group_by_at(selected_fill_parameter) %>% mutate(max_mean=max(mean)) %>% mutate(mean_norm=mean/max_mean) %>% mutate(max_sd=max(sd)) %>% mutate(sd_norm=sd/max_sd)
    topology_results_selected_analytical_norm_df = topology_results_selected_analytical_df %>% group_by_at(selected_fill_parameter) %>% mutate(model_result = model_result/max_value_in_dataset)
    predicted_results_norm_df = predicted_results_df %>% group_by_at(selected_fill_parameter) %>% mutate(model_result = model_result/max_value_in_dataset)
  } else if(topology_type_normalization == "divide.L"){
    if (is.null(selected_fill_parameter)){
      results_selected_norm_df = results_selected_df %>% inner_join((topology_results_selected_analytical_df %>% filter(model %in% topology_type_analytical_model) %>% select("L", "max_value_in_dataset", "type_dataset") %>% unique()), by="type_dataset") %>% mutate(norm = (get(boxplot_parameter)/max_value_in_dataset)/L) %>% rename(unnorm = boxplot_parameter) %>% rename(!!boxplot_parameter := norm)
      topology_results_selected_by_size_norm_df = topology_results_selected_by_size_df %>% inner_join((topology_results_selected_analytical_df %>% filter(model %in% topology_type_analytical_model) %>% select("L", "max_value_in_dataset", "type_dataset") %>% unique()), by="type_dataset") %>% group_by_at(selected_fill_parameter) %>% mutate(mean_norm = (mean/max_value_in_dataset)/L)
      predicted_results_norm_df = predicted_results_df %>% inner_join((topology_results_selected_analytical_df %>% filter(model %in% topology_type_analytical_model) %>% select("L", "type_dataset") %>% unique()), by="type_dataset") %>% group_by_at(selected_fill_parameter) %>% mutate(model_result = (model_result / max_value_in_dataset)/L)
      topology_results_selected_analytical_norm_df = topology_results_selected_analytical_df %>% group_by_at(selected_fill_parameter) %>% mutate(norm = (model_result/max_value_in_dataset)/L) %>% rename(unnorm = model_result) %>% rename(model_result = norm)
    } else{
      results_selected_norm_df = results_selected_df %>% inner_join((topology_results_selected_analytical_df %>% filter(model %in% topology_type_analytical_model) %>% select("L", "max_value_in_dataset", !!selected_fill_parameter) %>% unique()), by=selected_fill_parameter) %>% group_by_at(selected_fill_parameter) %>% mutate(norm = (get(boxplot_parameter)/max_value_in_dataset)/L) %>% rename(unnorm = boxplot_parameter) %>% rename(!!boxplot_parameter := norm)
      topology_results_selected_by_size_norm_df = topology_results_selected_by_size_df %>% inner_join((topology_results_selected_analytical_df %>% filter(model %in% topology_type_analytical_model) %>% select("L", "max_value_in_dataset", !!selected_fill_parameter) %>% unique()), by=selected_fill_parameter) %>% group_by_at(selected_fill_parameter) %>% mutate(mean_norm = (mean/max_value_in_dataset)/L)
      predicted_results_norm_df = predicted_results_df %>% inner_join((topology_results_selected_analytical_df %>% filter(model %in% topology_type_analytical_model) %>% select("L", !!selected_fill_parameter) %>% unique()), by=selected_fill_parameter) %>% group_by_at(selected_fill_parameter) %>% mutate(model_result = (model_result / max_value_in_dataset)/L)
      topology_results_selected_analytical_norm_df = topology_results_selected_analytical_df %>% group_by_at(selected_fill_parameter) %>% mutate(norm = (model_result/max_value_in_dataset)/L) %>% rename(unnorm = model_result) %>% rename(model_result = norm)
    }
  } else if(topology_type_normalization == "divide.max.num.links"){
    results_selected_norm_df = results_selected_df %>% 
      #inner_join(numbers_complete_graph_df, by = "type_dataset") %>% 
      group_by_at(selected_fill_parameter) %>% mutate(norm = get(boxplot_parameter)/total_num_edges) %>% rename(unnorm = boxplot_parameter) %>% rename(!!boxplot_parameter := norm)
    topology_results_selected_by_size_norm_df = topology_results_selected_by_size_df %>% 
      #inner_join(numbers_complete_graph_df, by = "type_dataset") %>% 
      group_by_at(selected_fill_parameter) %>% mutate(mean_norm=mean/total_num_edges)
    topology_results_selected_analytical_norm_df = topology_results_selected_analytical_df %>% inner_join(numbers_complete_graph_df, by = "type_dataset") %>% group_by_at(selected_fill_parameter) %>% mutate(norm = model_result/total_num_edges) %>% rename(unnorm = model_result) %>% rename(model_result = norm)
    predicted_results_norm_df = predicted_results_df %>% inner_join(numbers_complete_graph_df, by = "type_dataset") %>% group_by_at(selected_fill_parameter) %>% mutate(norm = model_result/total_num_edges) %>% rename(unnorm = model_result) %>% rename(model_result = norm)
    # print size for percentage of edges
    for(type_dataset_selected in unique(predicted_results_norm_df$type_dataset)){
      predicted_results_specific = predicted_results_norm_df %>% filter((model %in% c(topology_type_analytical_model)) & (type_dataset==type_dataset_selected))
      #print(predicted_results_specific[which.min(abs(0.25-predicted_results_specific$model_result)),])
      #print(predicted_results_specific[which.min(abs(0.50-predicted_results_specific$model_result)),])
      #print(predicted_results_specific[which.min(abs(0.75-predicted_results_specific$model_result)),])
    }
  }
  combination_metric = paste(combination_metric, "norm", sep="_")
}


#--------------------------#
# Write normalized results #
#--------------------------#

topology_results_file = paste(output_dir, '/', 'topology_results_norm_', topology_type_normalization, '_', type_data_selection, '.txt', sep='')
results_selected_norm_df %>% fwrite(topology_results_file)
topology_results_by_size_file = paste(output_dir, '/', 'topology_results_mean_norm_', topology_type_normalization, '_', type_data_selection, '.txt', sep='')
topology_results_selected_by_size_norm_df %>% fwrite(topology_results_by_size_file)
predictions_file = paste(output_dir, '/', 'predictions_', topology_type_normalization, '_', type_data_selection, '.txt', sep='')
predicted_results_norm_df %>% fwrite(predictions_file)
analytical_results_file = paste(output_dir, '/', 'analytical_model_results_', topology_type_normalization, '_', type_data_selection, '.txt', sep='')
topology_results_selected_analytical_norm_df %>% fwrite(analytical_results_file)

