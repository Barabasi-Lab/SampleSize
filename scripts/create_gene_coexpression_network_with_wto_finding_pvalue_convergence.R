#!/usr/bin/env Rscript
library(optparse)
require(data.table)
require(dplyr)
require(magrittr)
require(tidyr)
require(wTO)
set.seed(1510)

#### READ ARGUMENTS ####
option_list = list(
  make_option(c("-s", "--samples_file"), type="character", 
              help="Samples list file. A file with the samples that will be used from the RNAseq file, each one in a different line.", metavar="character"),
  make_option(c("-f", "--rnaseq_file"), type="character", 
              help="RNAseq file name", metavar="character"),
  make_option(c("-o", "--output_wto_network"), type="character", default="gene_coexpression_network_wto.net", 
              help="output file [default= %default]", metavar="character"),
  make_option(c("-b", "--output_wto_bootstrap"), type="character", default="gene_coexpression_network_bootstrap.csv", 
              help="output file [default= %default]", metavar="character"),
  make_option(c("-u", "--output_corr_network"), type="character", default="gene_coexpression_network_pearson.net", 
              help="output file [default= %default]", metavar="character"),
  make_option(c("-n", "--n_max_iterations"), type="integer", default=300, 
              help="Number of maximum wTO bootstrap iterations to calculate the p-value [default= %default]", metavar="integer"),
  make_option(c("-e", "--n_edges_to_check"), type="integer", default=50, 
              help="Number of edges to check for convergence of the p-value [default= %default]", metavar="integer"),
  make_option(c("-d", "--delta"), type="double", default=0.2,
              help="Value that defines the interval of confidence from which the p-values of the bootstrap repetitions are calculated [default= %default]", metavar="double"),
  make_option(c("-p", "--n_pvals_to_check"), type="integer", default=10,
              help="Number of p-values to check if they converge [default= %default]", metavar="integer"),
  make_option(c("-c", "--pval_cutoff"), type="double", default=0.01,
              help="Cut-off of the standard deviation of the p-values checked. If the standard deviation of the p-values is below the cut-off, it means we have reached convergence [default= %default]", metavar="double")
); 
# Example of execution
# Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/create_gene_coexpression_network_with_wto_finding_pvalue_convergence.R -s /home/j.aguirreplans/Projects/Scipher/SampleSize/data/sampling/GTEx/sampling_with_repetition/Whole.Blood_female/RNAseq_samples_Whole.Blood_female_size_10_rep_1.txt -f /home/j.aguirreplans/Databases/GTEx/v8/tpm_filtered_files_by_tissue/GTEx_RNASeq_Whole.Blood_female.gct -o /scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/test/wto_Whole.Blood_female_size_10_rep_1.net -b /scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/test/wto_boot_Whole.Blood_female_size_10_rep_1.net -u /scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/test/pearson_Whole.Blood_female_size_10_rep_1.net -n 300 -e 50 -d 0.2 -p 10 -c 0.01

# Read arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Check for missing arguments
if (is.null(opt$rnaseq_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (rnaseq file).n", call.=FALSE)
}
samples_file = opt$samples_file
rnaseq_file = opt$rnaseq_file
output_wto_network = opt$output_wto_network
output_wto_bootstrap = opt$output_wto_bootstrap
output_corr_network = opt$output_corr_network
n_max_iterations = strtoi(opt$n_max_iterations)
n_edges_to_check = strtoi(opt$n_edges_to_check)
delta = as.double(opt$delta)
n_pvals_to_check = strtoi(opt$n_pvals_to_check)
pval_cutoff = as.double(opt$pval_cutoff)

#samples_file = '/home/j.aguirreplans/Projects/Scipher/SampleSize/data/sampling/GTEx/sampling_with_repetition/Whole.Blood_female/RNAseq_samples_Whole.Blood_female_size_10_rep_1.txt'
#rnaseq_file = '/home/j.aguirreplans/Databases/GTEx/v8/tpm_filtered_files_by_tissue/GTEx_RNASeq_Whole.Blood_female.gct'
#output_wto_network = '/scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/test/wto_Whole.Blood_female_size_10_rep_1.net'
#output_wto_bootstrap = '/scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/test/wto_boot_Whole.Blood_female_size_10_rep_1.csv'
#output_corr_network = '/scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/test/pearson_Whole.Blood_female_size_10_rep_1.net'
#n_max_iterations = 300
#n_edges_to_check = 50
#delta = 0.2
#n_pvals_to_check = 10
#pval_cutoff = 0.01

# Get scripts dir
#initial.options <- commandArgs(trailingOnly = FALSE)
#scripts.dir <- dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))


#### READ DATASET ####

# Read samples
subsample = fread(samples_file)[, 1][[1]]

# Read RNAseq dataset (rows = genes, columns = samples)
rnaseq = fread(rnaseq_file)

# Subset gene expression datast by samples in the samples file
rnaseq = rnaseq %>% select(c(colnames(rnaseq)[1], all_of(subsample)))

# Remove the samples column and include it as "rownames"
gene.ids <- rnaseq[, 1][[1]]
rnaseq <- as.matrix(rnaseq[, -c(1)])
rownames(rnaseq) <- gene.ids

#rnaseq = rnaseq[row.names(rnaseq) %in% sample(row.names(rnaseq), size=2000, replace=FALSE),] # Check example with less genes


#### RUN CO-EXPRESSION WITH WTO CHECKING P-VALUE CONVERGENCE IN EACH ITERATION ####

# Define WTO function
wTO = function(A){
  A_TF = A
  C = as.matrix(A) %*% t(A)
  W = C + A_TF
  K = matrix(NA, nrow(A_TF), ncol(A_TF))
  KI = rowSums(abs(A), na.rm = T)
  for (ii in 1:nrow(A_TF)) {
    for (jj in 1:ncol(A_TF)) {
      K[ii, jj] = min(KI[ii], KI[jj])
    }
  }
  WTO = round(W/(K + 1 - abs(A_TF)), 3)
  return(WTO)
}


# Define function to calculate if a score in a bootstrap iteration is close to the real score
calculate_if_score_falls_into_interval = function(wTO.col, Iteration.wTO.col, delta){
  return(ifelse( Iteration.wTO.col > (wTO.col - delta) & Iteration.wTO.col < (wTO.col + delta), 0, 1))
}


# Define function to check convergence for multiple edges
check_pvalue_convergence_for_multiple_edges = function(wto_results_df, edges_for_convergence_checking, delta=0.2, pval_cutoff=0.01, n_pvals_to_check=5){
  
  edges_with_sd_below_cutoff = c()
  edges_with_diff_below_cutoff = c()
  if ((ncol(wto_results_df)-3) >= n_pvals_to_check){
    for (edge in edges_for_convergence_checking){
      
      # Get the expression of selected edge for all bootstrap repetitions
      edge_results = wto_results_df[, -c(1:3)] %>% filter(row.names(wto_results_df) == edge) %>% as.numeric() %>% data.frame()
      names(edge_results) = "wTO"
      edge_results$run = 1:nrow(edge_results)
      
      # To calculate the p-value, we add and substract the value of delta to the edge 1 wTO score without bootstrapping (C[1,3]). 
      edge_wto = (wto_results_df %>% filter(row.names(wto_results_df) == edge))$wTO
      min = (edge_wto - delta) %>% as.numeric()
      max = (edge_wto + delta) %>% as.numeric()
      
      # If the value is outside the interval obtained as a result of adding/substracting delta, it puts a 1
      edge_results$p = ifelse( edge_results$wTO > min & edge_results$wTO < max, 0, 1)
      
      # Here, we sum all the scores outside the limits defined by delta, to calculate the p-value
      edge_results$prop = cumsum(edge_results$p)
      edge_results$pval = edge_results$prop/edge_results$run
      
      # We get the last n p-values to check
      pvals_to_check = edge_results$pval[(nrow(edge_results)-n_pvals_to_check+1):nrow(edge_results)]
      
      # We calculate the sd between the last n p-values
      sd_pvals = sd(pvals_to_check)
      
      # We check if the sd is below the cut-off established
      pval_sd_below_cutoff = sd_pvals <= pval_cutoff
      edges_with_sd_below_cutoff = c(edges_with_sd_below_cutoff, pval_sd_below_cutoff)
      
      # We calculate the difference between consecutive pairs of p-values in the last n p-values
      diff_pvals = diff(pvals_to_check)

      # We check if the sd is below the cut-off established
      pval_diff_below_cutoff = all(diff_pvals <= pval_cutoff)
      edges_with_diff_below_cutoff = c(edges_with_diff_below_cutoff, pval_diff_below_cutoff)
      
    }
    return(list("edges_with_sd_below_cutoff" = edges_with_sd_below_cutoff, "edges_with_diff_below_cutoff" = edges_with_diff_below_cutoff))
  } else {
    return(list("edges_with_sd_below_cutoff" = c(FALSE), "edges_with_diff_below_cutoff" = c(FALSE)))
  }
  
}


# Define the function to plot the convergence of an edge
plot_pvalue_convergence_for_edge = function(wto_results_df, edge_rowname, delta=0.2){
  
  # Get the expression of selected edge for all bootstrap repetitions
  edge_results = wto_results_df[, -c(1:3)] %>% filter(row.names(wto_results_df) == edge_rowname) %>% as.numeric() %>% data.frame()
  names(edge_results) = "wTO"
  edge_results$run = 1:nrow(edge_results)
  
  # To calculate the p-value, we add and substract the value of delta to the edge 1 wTO score without bootstrapping (C[1,3]). 
  edge_wto = (wto_results_df %>% filter(row.names(wto_results_df) == edge_rowname))$wTO
  min = (edge_wto - delta) %>% as.numeric()
  max = (edge_wto + delta) %>% as.numeric()
  
  # If the value is outside the interval obtained as a result of adding/substracting delta, it puts a 1
  edge_results$p = ifelse( edge_results$wTO > min & edge_results$wTO < max, 0, 1)
  
  # Here, we sum all the scores outside the limits defined by delta, to calculate the p-value
  edge_results$prop = cumsum(edge_results$p)
  edge_results$pval = edge_results$prop/edge_results$run
  
  # Plot the edge p-values
  ggplot(edge_results) +
    aes(x = run, y = pval) +
    geom_point(size = 1L, colour = "#bd3786") +
    theme_minimal()
  
}


# Calculate Pearson correlation
A = CorrelationOverlap(Data = rnaseq, Overlap = row.names(rnaseq),
                       method = "p") %>% as.data.frame() 
A_net = A %>% wTO.in.line() %>% rename(Pearson=wTO)

# Save correlation in a file
fwrite(A_net, output_corr_network)
rm(A_net)

# Calculate wTO scores for the samples
wto_results_df = wTO(A) %>% wTO.in.line()

# Get the edges that we will check for convergence during iterations
edges_for_convergence_checking = sample(row.names(wto_results_df), size=n_edges_to_check, replace=FALSE)

# Start bootstrap iterations
for (i in 1:n_max_iterations){
  
  print(i)
  A = CorrelationOverlap(Data =  rnaseq[, sample(1:ncol(rnaseq),
                                                      replace = TRUE)], Overlap = row.names(rnaseq),
                              method = "p") %>% as.data.frame()
  wto_results_df <- inner_join(wto_results_df, wTO(A) %>% wTO.in.line() %>% rename_with(.fn = ~paste0("C_", i), .cols = "wTO"), by=c("Node.1" = "Node.1", "Node.2" = "Node.2"))
  results_below_butoff = check_pvalue_convergence_for_multiple_edges(wto_results_df, edges_for_convergence_checking, delta, pval_cutoff, n_pvals_to_check)
  #print(results_below_butoff$edges_with_sd_below_cutoff)
  #print(length(results_below_butoff$edges_with_sd_below_cutoff[results_below_butoff$edges_with_sd_below_cutoff == TRUE]))
  print(results_below_butoff$edges_with_diff_below_cutoff)
  print(length(results_below_butoff$edges_with_diff_below_cutoff[results_below_butoff$edges_with_diff_below_cutoff == TRUE]))
  if (all(results_below_butoff$edges_with_diff_below_cutoff)){
    break
  }
  
}

# Write the bootstrap result
fwrite(wto_results_df, output_wto_bootstrap)

rm(rnaseq)
rm(A)


# Calculate the p-value of the network
wto_network = wto_results_df %>%
  pivot_longer(
    cols = starts_with("C_"),
    names_to = "Iteration",
    names_prefix = "C_",
    values_to = "Iteration.wTO"
  ) %>% 
  mutate(across(Iteration, as.numeric)) %>% 
  group_by(Node.1, Node.2, wTO) %>% 
  mutate(p = calculate_if_score_falls_into_interval(wTO, Iteration.wTO, delta), prop = cumsum(p), pval = prop/Iteration) %>% 
  #mutate(p = calculate_if_score_falls_into_interval(wTO, Iteration.wTO, delta), prop = cumsum(p), pval = prop/Iteration, padj = p.adjust(pval, method = "BH")) %>% 
  filter(Iteration == max(Iteration)) %>%
  select(Node.1, Node.2, wTO, pval) %>%
  #select(Node.1, Node.2, wTO, pval, padj) %>%
  ungroup()
wto_network$padj = p.adjust(wto_network$pval, method = "BH")

rm(wto_results_df)

# Write the final network
fwrite(wto_network, output_wto_network)

rm(wto_network)








#plot_pvalue_convergence_for_edge(wto_results_df, edge_rowname=edges_for_convergence_checking[1], delta=0.2)
#plot_pvalue_convergence_for_edge(wto_results_df, edge_rowname=edges_for_convergence_checking[2], delta=0.2)

# Define function to check convergence for multiple edges
# I could not finish it because I don't find a way to calculate the difference...
check_pvalue_convergence_for_multiple_edges_without_loop = function(wto_results_df, edges_for_convergence_checking, delta=0.2, pval_cutoff=0.01, n_pvals_to_check=5){
  
  edges_with_sd_below_cutoff = c()
  edges_with_diff_below_cutoff = c()
  if ((ncol(wto_results_df)-3) >= n_pvals_to_check){
    
    wto_pval_results = wto_results_df %>%
      filter(row.names(wto_results_df) %in% edges_for_convergence_checking) %>%
      pivot_longer(
        cols = starts_with("C_"),
        names_to = "Iteration",
        names_prefix = "C_",
        values_to = "Iteration.wTO"
      ) %>% 
      # Transform iteration to numeric
      mutate(across(Iteration, as.numeric)) %>% 
      # Group by edge
      group_by(Node.1, Node.2, wTO) %>% 
      # Calculate p-value
      mutate(p = calculate_if_score_falls_into_interval(wTO, Iteration.wTO, delta), prop = cumsum(p), pval = prop/Iteration) %>% 
      # Get the last n p-values for each edge
      filter(Iteration %in% ((max(Iteration)-n_pvals_to_check+1):max(Iteration))) %>%
      mutate(sd = sd(pval)) %>%
      select(Node.1, Node.2, wTO, Iteration, pval) %>%
      ungroup()

    return(list("edges_with_sd_below_cutoff" = edges_with_sd_below_cutoff, "edges_with_diff_below_cutoff" = edges_with_diff_below_cutoff))
  } else {
    return(list("edges_with_sd_below_cutoff" = c(FALSE), "edges_with_diff_below_cutoff" = c(FALSE)))
  }

}


