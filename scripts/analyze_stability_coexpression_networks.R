#!/usr/bin/env Rscript
packrat::init("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeR")
library(optparse)
library(data.table)
library(dplyr)
library(ggplot2)
library(igraph)
library(mltools)
library(wTO)
set.seed(1510)

#### READ ARGUMENTS ####
option_list = list(
  make_option(c("-c", "--coexpression_network_file"), type="character", 
              help="Co-expression network from limited samples file", metavar="character"),
  make_option(c("-a", "--coexpression_network_all_samples_file"), type="character", 
              help="Co-expression network from all samples file", metavar="character"),
  make_option(c("-p", "--ppi_file"), type="character", 
              help="Protein-protein interaction network file.", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default="analyze_stability_coexpression_networks.txt", 
              help="Output file [default= %default]", metavar="character")
  #make_option(c("-o", "--output_dir"), type="character", default="analyze_stability_coexpression_networks.txt", 
  #            help="Output dir [default= %default]", metavar="character")
); 
# Example of execution
# Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_stability_coexpression_networks.R -c /home/j.aguirreplans/data/PPI/interactome_tissue_specific/interactome_2019_Spleen_female.csv -p /home/j.aguirreplans/data/PPI/interactome_tissue_specific/interactome_2019_Spleen_female.csv -o /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_GTEx/Spleen_female/wgcna_RNAseq_samples_Spleen_female_all_samples.net/comparison_coexpression_ppi.txt

# Read arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Check for missing arguments
if ((is.null(opt$coexpression_network_file)) | (is.null(opt$coexpression_network_all_samples_file)) | (is.null(opt$ppi_file))){
  print_help(opt_parser)
  stop("At least one argument must be supplied (coexpression_network_file, ppi_file).n", call.=FALSE)
}
tissue_sex_network_file = opt$coexpression_network_file
tissue_sex_network_all_samples_file = opt$coexpression_network_all_samples_file
ppi_tissue_sex_filt_file = opt$ppi_file
#output_dir = opt$output_dir
output_file = opt$output_file

#tissue_sex_network_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/Spleen_female/wgcna_RNAseq_samples_Spleen_female_size_60_rep_8.net"
#tissue_sex_network_all_samples_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/Spleen_female/wgcna_RNAseq_samples_Spleen_female_all_samples.net"
#ppi_tissue_sex_filt_file = "/home/j.aguirreplans/data/PPI/interactome_tissue_specific/interactome_2019_Spleen_female.csv"
#output_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_GTEx/Spleen_female/wto_RNAseq_samples_Spleen_female_all_samples.net"
#output_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_GTEx/Spleen_female/comparison_coexpression_ppi_wto_RNAseq_samples_Spleen_female_all_samples.net"

# Read gene co-expression networks
read_coexpression_network <- function(coexpression_network_file){
  coexpression_network_df = as.data.frame(fread(coexpression_network_file))
  if (!(grepl('wto', coexpression_network_file, fixed = TRUE))){
    # Pass the matrix to the wTO in-line format
    rownames(coexpression_network_df) = colnames(coexpression_network_df)
    coexpression_network_df <- coexpression_network_df[order(rownames(coexpression_network_df)), order(colnames(coexpression_network_df))]
    coexpression_network_df = coexpression_network_df %>% wTO.in.line() %>% rename(score=wTO)
  } else {
    coexpression_network_df = coexpression_network_df %>% rename(score=wTO)
    coexpression_network_df = coexpression_network_df %>% filter(pval.adj < 0.05)
  }
}

tissue_sex_network_df = read_coexpression_network(tissue_sex_network_file) %>% rename("score.subset.network" = "score")
tissue_sex_network_all_samples_df = read_coexpression_network(tissue_sex_network_all_samples_file) %>% rename("score.all.samples.network" = "score")

# Merge two dataframes
tissue_sex_network_merged_df <- inner_join(tissue_sex_network_df, tissue_sex_network_all_samples_df, by=c("Node.1" = "Node.1", "Node.2" = "Node.2"))
rm(tissue_sex_network_df)
rm(tissue_sex_network_all_samples_df)


# Read PPI network
ppi_tissue_sex_filt_df = fread(ppi_tissue_sex_filt_file) %>% select(proteinA_symbol, proteinB_symbol)
ppi_tissue_sex_filt_net = graph_from_data_frame(ppi_tissue_sex_filt_df, directed=F) %>% simplify()

# Calculate clustering coefficient, degree and betweenness centrality of nodes
# Keep it in a dataframe
node_topology_info_df = data.frame(nodes=V(ppi_tissue_sex_filt_net)$name, 
                                   local_clustering_coefficient=transitivity(ppi_tissue_sex_filt_net, type="local", isolate="zero"), 
                                   degree_centrality = as.vector(degree(ppi_tissue_sex_filt_net)),
                                   betweenness_centrality = as.vector(betweenness(ppi_tissue_sex_filt_net, directed=FALSE))
)

# Calculate distances between pairs of nodes of the PPI network
ppi_tissue_sex_filt_distances = distances(ppi_tissue_sex_filt_net)
ppi_tissue_sex_filt_distances <- ppi_tissue_sex_filt_distances[order(rownames(ppi_tissue_sex_filt_distances)), order(colnames(ppi_tissue_sex_filt_distances))]
# Convert the matrix to an in-line format and remove infinite distances (from different components)
ppi_tissue_sex_filt_distances = ppi_tissue_sex_filt_distances %>% wTO.in.line() %>% rename(distances=wTO)
ppi_tissue_sex_filt_distances = ppi_tissue_sex_filt_distances %>% filter(!(is.infinite(distances)))
# Merge node topology properties with distances
ppi_tissue_sex_filt_distances = left_join(ppi_tissue_sex_filt_distances, node_topology_info_df, by=c("Node.1" = "nodes")) %>% rename(local_clustering_coefficient.1=local_clustering_coefficient, degree_centrality.1=degree_centrality, betweenness_centrality.1=betweenness_centrality)
ppi_tissue_sex_filt_distances = left_join(ppi_tissue_sex_filt_distances, node_topology_info_df, by=c("Node.2" = "nodes")) %>% rename(local_clustering_coefficient.2=local_clustering_coefficient, degree_centrality.2=degree_centrality, betweenness_centrality.2=betweenness_centrality)
rm(node_topology_info_df)
ppi_tissue_sex_filt_distances$local_clustering_coefficient.mean = rowMeans((ppi_tissue_sex_filt_distances %>% select(local_clustering_coefficient.1, local_clustering_coefficient.2)))
ppi_tissue_sex_filt_distances$local_clustering_coefficient.1 = NULL
ppi_tissue_sex_filt_distances$local_clustering_coefficient.2 = NULL
ppi_tissue_sex_filt_distances$degree_centrality.mean = rowMeans((ppi_tissue_sex_filt_distances %>% select(degree_centrality.1, degree_centrality.2)))
ppi_tissue_sex_filt_distances$degree_centrality.1 = NULL
ppi_tissue_sex_filt_distances$degree_centrality.2 = NULL
ppi_tissue_sex_filt_distances$betweenness_centrality.mean = rowMeans((ppi_tissue_sex_filt_distances %>% select(betweenness_centrality.1, betweenness_centrality.2)))
ppi_tissue_sex_filt_distances$betweenness_centrality.1 = NULL
ppi_tissue_sex_filt_distances$betweenness_centrality.2 = NULL

# Merge the two dataframes (gene co-expression network and distances network) as both have all possible combinations of pairs of nodes
ppi_gene_coexpression_merged_df <- inner_join(ppi_tissue_sex_filt_distances, tissue_sex_network_merged_df, by=c("Node.1" = "Node.1", "Node.2" = "Node.2"))
rm(ppi_tissue_sex_filt_distances)
rm(tissue_sex_network_merged_df)

# Calculate difference between networks
ppi_gene_coexpression_merged_df$difference = ppi_gene_coexpression_merged_df$score.all.samples.network - ppi_gene_coexpression_merged_df$score.subset.network

# Get parameters from different ranges of differences betweens cores
difference_df = data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("left.range", "right.range", "number.of.edges"))))
ranges_df = data.frame(left.range=c(-Inf, -0.1, -0.05, -0.01, 0, 0.01, 0.05, 0.1), right.range=c(-0.1, -0.05, -0.01, 0, 0.01, 0.05, 0.1, Inf))
for(range.id in rownames(ranges_df)){
  left.range = ranges_df[range.id,"left.range"]
  right.range = ranges_df[range.id,"right.range"]
  subset_df = ppi_gene_coexpression_merged_df %>% filter((difference >= left.range) & (difference < right.range))
  difference_df = rbind(difference_df, data.frame(left.range=left.range, right.range=right.range, 
                                                  number.of.edges=length(subset_df$difference),
                                                  difference.mean=mean(subset_df$difference),
                                                  score.subset.mean=mean(subset_df$score.subset.network),
                                                  score.all.samples.mean=mean(subset_df$score.all.samples.network),
                                                  distance.mean=mean(subset_df$distances),
                                                  local_clustering_coefficient.mean=mean(subset_df$local_clustering_coefficient.mean),
                                                  degree_centrality.mean=mean(subset_df$degree_centrality.mean),
                                                  betweenness_centrality.mean=mean(subset_df$betweenness_centrality.mean)
                                                  ))
}
rm(subset_df)

scores_df = data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("left.range", "right.range", "number.of.edges"))))
scores_direct_df = data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("left.range", "right.range", "number.of.edges"))))
ranges_df = data.frame(left.range=seq(-1, 0.9, 0.1), right.range=seq(-0.9, 1, 0.1))
for(range.id in rownames(ranges_df)){
  left.range = ranges_df[range.id,"left.range"]
  right.range = ranges_df[range.id,"right.range"]
  subset_df = ppi_gene_coexpression_merged_df %>% filter((score.subset.network >= left.range) & (score.subset.network < right.range))
  scores_df = rbind(scores_df, data.frame(left.range=left.range, right.range=right.range, 
                                          number.of.edges=length(subset_df$difference),
                                          difference.mean=mean(subset_df$difference),
                                          score.subset.mean=mean(subset_df$score.subset.network),
                                          score.all.samples.mean=mean(subset_df$score.all.samples.network),
                                          distance.mean=mean(subset_df$distances),
                                          local_clustering_coefficient.mean=mean(subset_df$local_clustering_coefficient.mean),
                                          degree_centrality.mean=mean(subset_df$degree_centrality.mean),
                                          betweenness_centrality.mean=mean(subset_df$betweenness_centrality.mean)
  ))
  subset_df = ppi_gene_coexpression_merged_df %>% filter((distances < 2) & (score.subset.network >= left.range) & (score.subset.network < right.range))
  scores_direct_df = rbind(scores_direct_df, data.frame(left.range=left.range, right.range=right.range, 
                                                        number.of.edges=length(subset_df$difference),
                                                        difference.mean=mean(subset_df$difference),
                                                        score.subset.mean=mean(subset_df$score.subset.network),
                                                        score.all.samples.mean=mean(subset_df$score.all.samples.network),
                                                        distance.mean=mean(subset_df$distances),
                                                        local_clustering_coefficient.mean=mean(subset_df$local_clustering_coefficient.mean),
                                                        degree_centrality.mean=mean(subset_df$degree_centrality.mean),
                                                        betweenness_centrality.mean=mean(subset_df$betweenness_centrality.mean)
  ))
}
rm(subset_df)

# Create output file
difference_df %>% fwrite(output_file)




