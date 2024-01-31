#### READ ARGUMENTS ####
option_list = list(
  make_option(c("-c", "--coexpression_network_file"), type="character", 
              help="Co-expression network file", metavar="character"),
  make_option(c("-a", "--disease_genes_file"), type="character", default = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/disease_genes/disease_genes_info.csv",
              help="File containing the disease genes in the expression dataset", metavar="character"),
  make_option(c("-b", "--essential_genes_file"), type="character", default = "/home/j.aguirreplans/Databases/OGEE/OGEE_esential_genes.csv",
              help="File containing the essential genes in the expression dataset", metavar="character"),
  #make_option(c("-n", "--output_file_name"), type="character", default="network", 
  #            help="Output file name [default= %default]", metavar="character"),
  #make_option(c("-o", "--output_dir"), type="character", default="/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_scipher/all_samples", 
  #            help="Output dir [default= %default]", metavar="character")
  make_option(c("-t", "--output_topology_file"), type="character", default="analysis_topology.txt", 
              help="Output topology file [default= %default]", metavar="character"),
  make_option(c("-d", "--output_disease_genes_file"), type="character", default="analysis_essential_genes.txt", 
              help="Output disease genes file [default= %default]", metavar="character"),
  make_option(c("-e", "--output_essential_genes_file"), type="character", default="analysis_disease_genes.txt", 
              help="Output file [default= %default]", metavar="character")
); 
# Example of execution
# Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_coexpression_network_by_top_scoring_edges.R -c /home/j.aguirreplans/data/PPI/interactome_tissue_specific/interactome_2019_Spleen_female.csv 

# Read arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Check for missing arguments
if (is.null(opt$coexpression_network_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (coexpression_network_file, ppi_file).n", call.=FALSE)
}
coexpression_network_file = opt$coexpression_network_file
disease_genes_file = opt$disease_genes_file
essential_genes_file = opt$essential_genes_file
output_file_name = opt$output_file_name
output_topology_file = opt$output_topology_file
output_essential_genes_file = opt$output_essential_genes_file
output_disease_genes_file = opt$output_disease_genes_file

#coexpression_network_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_scipher/all_samples/wgcna_scipher_all_samples_size_100_rep_1.net"
#output_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_scipher/all_samples/analysis_by_top_scoring_edges_wgcna_scipher_all_samples_size_100_rep_1.txt"
#disease_genes_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/disease_genes/disease_genes_info.csv"
#essential_genes_file = "/home/j.aguirreplans/Databases/OGEE/OGEE_esential_genes.csv"

# Read co-expression network
coexpression_df = as.data.frame(fread(coexpression_network_file))
if(grepl('wto', coexpression_network_file, fixed = TRUE)){
  score_threshold_list = c(0.3, 0.5, 0.7)
  coexpression_df = coexpression_df %>% rename("score"="wTO")
  coexpression_df = coexpression_df %>% filter(pval.adj < 0.05)
  #score_quantile = quantile(abs(ppi_gene_coexpression_merged_df$score), probs=c(0.8, 0.9))
} else if((grepl('spearman', coexpression_network_file, fixed = TRUE))){
  score_threshold_list = c(0.6, 0.65, 0.7, 0.75)
  coexpression_df = coexpression_df %>% rename("score"="spearman")
  coexpression_df = coexpression_df %>% filter(pvalue < 0.05)
} else if((grepl('pearson', coexpression_network_file, fixed = TRUE))){
  score_threshold_list = c(0.6, 0.65, 0.7, 0.75)
  coexpression_df = coexpression_df %>% rename("score"="pearson")
  coexpression_df = coexpression_df %>% filter(pvalue < 0.05)
} else if((grepl('wgcna', coexpression_network_file, fixed = TRUE))){
  require(WGCNA)
  score_threshold_list = c("hard", 0.3, 0.4, 0.5, 0.6)
  # Calculate hard threshold
  hard_threshold = pickHardThreshold(coexpression_df, RsquaredCut = 0.80)
  # Pass the matrix to the wTO in-line format
  rownames(coexpression_df) = colnames(coexpression_df)
  coexpression_df <- coexpression_df[order(rownames(coexpression_df)), order(colnames(coexpression_df))]
  coexpression_df = coexpression_df %>% wTO.in.line() %>% rename(score=wTO)
} else {
  score_threshold_list = c(0.3, 0.4, 0.5)
  # Pass the matrix to the wTO in-line format
  rownames(coexpression_df) = colnames(coexpression_df)
  coexpression_df <- coexpression_df[order(rownames(coexpression_df)), order(colnames(coexpression_df))]
  coexpression_df = coexpression_df %>% wTO.in.line() %>% rename(score=wTO)
}

genes_network = unique(c(as.vector(coexpression_df$Node.1), as.vector(coexpression_df$Node.2)))
disease_genes_df = fread(disease_genes_file)
disease_genes_df = disease_genes_df %>% filter(gene %in% genes_network)
essential_genes = unique(fread(essential_genes_file)$gene)
essential_genes = essential_genes[essential_genes %in% genes_network]

# Define result dataframes
topology_df = data.frame(matrix(ncol=14,nrow=0, dimnames=list(NULL, c("score_threshold", "num_nodes", "num_edges", "av_degree", "av_path_length", "av_clustering_coef", "num_components", "num_lcc_nodes", "num_lcc_edges", "lcc_z", "lcc_pvalue", "max_k", "num_main_core_nodes", "num_main_core_edges"))))
essential_gene_results_df = data.frame(matrix(ncol=9,nrow=0, dimnames=list(NULL, c("score_threshold", "num_essential_genes", "fraction_essential_genes", "num_components", "num_lcc_nodes", "fraction_essential_lcc_nodes", "num_lcc_edges", "lcc_z", "lcc_pvalue"))))
disease_gene_results_df <- setNames(data.frame(matrix(ncol = 11, nrow = 0)), c('score_threshold', 'disease', 'disease_class', 'num_disease_genes', 'fraction_disease_genes', 'num_disease_components', 'num_disease_lcc_nodes', 'fraction_disease_lcc_nodes', 'num_disease_lcc_edges', "disease_lcc_z", "disease_lcc_pvalue"))


for(score_threshold in score_threshold_list){
  
  # Select top scoring edges
  if(score_threshold == "hard"){
    if(!(is.na(hard_threshold$cutEstimate))){
      selected_coexpression_df = coexpression_df %>% filter(score > hard_threshold$cutEstimate)
    } else {
      next
    }
  } else {
    selected_coexpression_df = coexpression_df %>% filter(abs(score) > score_threshold)
  }
  selected_coexpression_net = graph_from_data_frame(selected_coexpression_df, directed=F) %>% simplify()
  
  # Calculate topological parameters
  num_nodes = gorder(selected_coexpression_net)
  num_edges = gsize(selected_coexpression_net)
  av_degree = mean(degree(selected_coexpression_net))
  av_path_length = mean_distance(selected_coexpression_net)
  av_clustering_coef = transitivity(selected_coexpression_net)
  
  # Components analysis
  components_net = components(selected_coexpression_net)
  lcc = induced.subgraph(selected_coexpression_net, vids = V(selected_coexpression_net)[components_net$membership == which.max(components_net$csize)] )
  num_lcc_nodes = gorder(lcc)
  num_lcc_edges = gsize(lcc)
  lcc_sig = LCC_Significance(N = 1000, 
                             Targets = V(lcc)$name, 
                             G = selected_coexpression_net)
  
  # K-core analysis
  k_core = coreness(selected_coexpression_net)
  max_k = max(k_core)
  main_core = induced.subgraph(selected_coexpression_net, vids = names(k_core[k_core==max_k]))
  num_main_core_nodes = gorder(main_core)
  num_main_core_edges = gsize(main_core)
  
  topology_df = rbind(topology_df, data.frame(score_threshold=score_threshold, num_nodes=num_nodes, num_edges=num_edges, av_degree=av_degree, av_path_length=av_path_length, av_clustering_coef=av_clustering_coef, num_components=components_net$no, num_lcc_nodes=num_lcc_nodes, num_lcc_edges=num_lcc_edges, lcc_z=lcc_sig$Z, lcc_pvalue=lcc_sig$emp_p, max_k=max_k, num_main_core_nodes=num_main_core_nodes, num_main_core_edges=num_main_core_edges))
  
  
  # Essential genes analysis
  essential_genes_in_network = essential_genes[essential_genes %in% V(selected_coexpression_net)$name]
  essential_subgraph = induced.subgraph(selected_coexpression_net, vids=essential_genes_in_network)
  essential_components = components(essential_subgraph)
  essential_lcc = induced.subgraph(essential_subgraph, vids = V(essential_subgraph)[essential_components$membership == which.max(essential_components$csize)] )
  num_essential_lcc_nodes = gorder(essential_lcc)
  num_essential_lcc_edges = gsize(essential_lcc)
  essential_lcc_sig = LCC_Significance(N = 1000, Targets = V(essential_lcc)$name, G = selected_coexpression_net)
  essential_gene_results_df = rbind(essential_gene_results_df, data.frame(score_threshold=score_threshold, num_essential_genes=length(essential_genes_in_network), fraction_essential_genes=length(essential_genes_in_network)/length(essential_genes), num_components=essential_components$no, num_essential_lcc_nodes=num_essential_lcc_nodes, fraction_essential_lcc_nodes=num_essential_lcc_nodes/length(essential_genes), num_lcc_edges=num_essential_lcc_edges, lcc_z=essential_lcc_sig$Z, lcc_pvalue=essential_lcc_sig$emp_p))
  
  # Disease genes analysis
  for (disease in sort(unique(disease_genes_df$disease))){
    disease_selected_df = disease_genes_df %>% filter(disease == !!disease)
    disease_selected_genes = unique(disease_selected_df$gene)
    disease_selected_genes_in_network = disease_selected_genes[disease_selected_genes %in% V(selected_coexpression_net)$name]
    if (length(disease_selected_genes) > 19){
      disease_subgraph = induced.subgraph(selected_coexpression_net, vids=disease_selected_genes_in_network)
      disease_components = components(disease_subgraph)
      disease_lcc = induced.subgraph(disease_subgraph, vids = V(disease_subgraph)[disease_components$membership == which.max(disease_components$csize)] )
      num_disease_lcc_nodes = gorder(disease_lcc)
      num_disease_lcc_edges = gsize(disease_lcc)
      disease_lcc_sig = LCC_Significance(N = 1000, Targets = V(disease_lcc)$name, G = selected_coexpression_net)
      disease_selected_classes = unique(disease_selected_df$root.concept)
      for (disease_selected_class in disease_selected_classes){
        disease_gene_results_df = rbind(disease_gene_results_df, data.frame(score_threshold=score_threshold, disease=disease, disease_class=disease_selected_class, num_disease_genes=length(disease_selected_genes_in_network), fraction_disease_genes=length(disease_selected_genes_in_network)/length(disease_selected_genes), num_disease_components=disease_components$no, num_disease_lcc_nodes=num_disease_lcc_nodes, fraction_disease_lcc_nodes=num_disease_lcc_nodes/length(disease_selected_genes), num_disease_lcc_edges=num_disease_lcc_edges, disease_lcc_z=disease_lcc_sig$Z, disease_lcc_pvalue=disease_lcc_sig$emp_p))
      }
    }
  }
  
  
}

# Create output file
topology_df %>% fwrite(output_topology_file)
essential_gene_results_df %>% fwrite(output_essential_genes_file)
disease_gene_results_df %>% fwrite(output_disease_genes_file)


