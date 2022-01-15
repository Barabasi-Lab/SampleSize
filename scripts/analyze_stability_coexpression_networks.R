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
  make_option(c("-d", "--difference_file"), type="character", default="difference.txt", 
              help="Output file [default= %default]", metavar="character"),
  make_option(c("-s", "--scores_file"), type="character", default="scores.txt", 
              help="Output file [default= %default]", metavar="character"),
  make_option(c("-t", "--threshold_results_file"), type="character", default="threshold_results.txt", 
              help="Output file [default= %default]", metavar="character"),
  make_option(c("-r", "--disease_results_file"), type="character", default="disease_results.txt", 
              help="Output file [default= %default]", metavar="character")
); 
# Example of execution
# Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_stability_coexpression_networks.R -c /scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/Spleen_female/wgcna_RNAseq_samples_Spleen_female_size_10_rep_2.net -a /scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/Spleen_female/wgcna_RNAseq_samples_Spleen_female_all_samples.net -p /home/j.aguirreplans/data/PPI/interactome_tissue_specific/interactome_2019_Spleen_female.csv -d /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_GTEx/Spleen_female/analysis_score_difference_wgcna_RNAseq_samples_Spleen_female_size_10_rep_2.txt -s /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_GTEx/Spleen_female/analysis_score_ranges_wgcna_RNAseq_samples_Spleen_female_size_10_rep_2.txt -t /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_GTEx/Spleen_female/analysis_score_thresholds_wgcna_RNAseq_samples_Spleen_female_size_10_rep_2.txt

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
difference_file = opt$difference_file
scores_file = opt$scores_file
threshold_results_file = opt$threshold_results_file
disease_results_file = opt$disease_results_file

#tissue_sex_network_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/Spleen_female/wgcna_RNAseq_samples_Spleen_female_size_60_rep_8.net"
#tissue_sex_network_all_samples_file = "/scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/Spleen_female/wgcna_RNAseq_samples_Spleen_female_all_samples.net"
#ppi_tissue_sex_filt_file = "/home/j.aguirreplans/data/PPI/interactome_tissue_specific/interactome_2019_Spleen_female.csv"

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



# Parse MESH
mesh_file = '/home/j.aguirreplans/Databases/MeSH/mtrees2021.bin'
mesh_df = fread(mesh_file, header=FALSE, col.names = c("concept", "concept.id"))
mesh_df$root.concept.id = substr(mesh_df$concept.id,1,3)
root.concepts = mesh_df %>% filter(concept.id==root.concept.id) %>% select(concept, concept.id) %>% rename("root.concept"="concept", "root.concept.id"="concept.id")
mesh_df = left_join(mesh_df, root.concepts, by=c("root.concept.id" = "root.concept.id"))
mesh_df$concept = tolower(mesh_df$concept)
mesh_df$root.concept = tolower(mesh_df$root.concept)

# Parse disease2gene
processDisease2Genes = function(filepath) {
  disease_genes_df=data.frame(matrix(ncol=2,nrow=0, dimnames=list(NULL, c("disease", "gene"))))
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    lineValues=strsplit(line, split="\t")[[1]]
    disease=lineValues[2]
    genes=lineValues[3:length(lineValues)]
    disease2genes=data.frame(disease=disease,gene=genes)
    disease_genes_df=rbind(disease_genes_df, disease2genes)
  }
  close(con)
  return(disease_genes_df)
}
disease_genes_file = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/Guney2016_GenesDisease.tsv"
disease_genes_df = processDisease2Genes(disease_genes_file)

# Merge disease2gene table with MESH
disease_genes_df = inner_join(disease_genes_df, mesh_df, by=c("disease" = "concept"))


# Read PPI network
ppi_tissue_sex_filt_df = fread(ppi_tissue_sex_filt_file)

# Merge disease2gene table with gene symbols
ppi_gene_info = rbind((ppi_tissue_sex_filt_df %>% select(proteinA, proteinA_symbol) %>% rename(gene=proteinA, gene_symbol=proteinA_symbol)), (ppi_tissue_sex_filt_df %>% select(proteinB, proteinB_symbol) %>% rename(gene=proteinB, gene_symbol=proteinB_symbol)))  %>% unique()
ppi_gene_info$gene = as.character(ppi_gene_info$gene)
disease_genes_df = inner_join(disease_genes_df, ppi_gene_info, by=c("gene" = "gene"))
rm(ppi_gene_info)

# Read network as igraph object
ppi_tissue_sex_filt_df = ppi_tissue_sex_filt_df %>% select(proteinA_symbol, proteinB_symbol)
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
ppi_tissue_sex_filt_distances = left_join(ppi_tissue_sex_filt_distances, node_topology_info_df, by=c("Node.1" = "nodes")) %>% rename("local_clustering_coefficient.1"="local_clustering_coefficient", "degree_centrality.1"="degree_centrality", "betweenness_centrality.1"="betweenness_centrality")
ppi_tissue_sex_filt_distances = left_join(ppi_tissue_sex_filt_distances, node_topology_info_df, by=c("Node.2" = "nodes")) %>% rename("local_clustering_coefficient.2"="local_clustering_coefficient", "degree_centrality.2"="degree_centrality", "betweenness_centrality.2"="betweenness_centrality")
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
ppi_gene_coexpression_merged_df$difference = abs(ppi_gene_coexpression_merged_df$score.all.samples.network - ppi_gene_coexpression_merged_df$score.subset.network)

# Get parameters from different ranges of differences betweens scores
difference_df = data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("left.range", "right.range", "number.of.edges", "type_interaction"))))
ranges_df = data.frame(left.range=c(0, 0.01, 0.05, 0.1, 0), right.range=c(0.01, 0.05, 0.1, Inf, Inf))
for(range.id in rownames(ranges_df)){
  left.range = ranges_df[range.id,"left.range"]
  right.range = ranges_df[range.id,"right.range"]
  subset_df = ppi_gene_coexpression_merged_df %>% filter((difference >= left.range) & (difference < right.range))
  subset_direct_df = subset_df %>% filter(distances < 2)
  subset_indirect_df = subset_df %>% filter(distances >= 2)
  difference_df = rbind(difference_df, data.frame(left.range=left.range, right.range=right.range, 
                                                  number.of.edges=length(subset_df$difference),
                                                  difference.mean=mean(subset_df$difference),
                                                  score.subset.mean=mean(subset_df$score.subset.network),
                                                  score.all.samples.mean=mean(subset_df$score.all.samples.network),
                                                  distance.mean=mean(subset_df$distances),
                                                  local_clustering_coefficient.mean=mean(subset_df$local_clustering_coefficient.mean),
                                                  degree_centrality.mean=mean(subset_df$degree_centrality.mean),
                                                  betweenness_centrality.mean=mean(subset_df$betweenness_centrality.mean),
                                                  type_interaction="all"
                                                  ))
  difference_df = rbind(difference_df, data.frame(left.range=left.range, right.range=right.range, 
                                                  number.of.edges=length(subset_direct_df$difference),
                                                  difference.mean=mean(subset_direct_df$difference),
                                                  score.subset.mean=mean(subset_direct_df$score.subset.network),
                                                  score.all.samples.mean=mean(subset_direct_df$score.all.samples.network),
                                                  distance.mean=mean(subset_direct_df$distances),
                                                  local_clustering_coefficient.mean=mean(subset_direct_df$local_clustering_coefficient.mean),
                                                  degree_centrality.mean=mean(subset_direct_df$degree_centrality.mean),
                                                  betweenness_centrality.mean=mean(subset_direct_df$betweenness_centrality.mean),
                                                  type_interaction="direct"
                                                  ))
  difference_df = rbind(difference_df, data.frame(left.range=left.range, right.range=right.range, 
                                                  number.of.edges=length(subset_indirect_df$difference),
                                                  difference.mean=mean(subset_indirect_df$difference),
                                                  score.subset.mean=mean(subset_indirect_df$score.subset.network),
                                                  score.all.samples.mean=mean(subset_indirect_df$score.all.samples.network),
                                                  distance.mean=mean(subset_indirect_df$distances),
                                                  local_clustering_coefficient.mean=mean(subset_indirect_df$local_clustering_coefficient.mean),
                                                  degree_centrality.mean=mean(subset_indirect_df$degree_centrality.mean),
                                                  betweenness_centrality.mean=mean(subset_indirect_df$betweenness_centrality.mean),
                                                  type_interaction="indirect"
                                                  ))
}


# Create output file
difference_df %>% fwrite(difference_file)


# Calculate parameters depending on range of scores
scores_df = data.frame(matrix(ncol=11,nrow=0, dimnames=list(NULL, c("left.range", "right.range", "number.of.edges", "difference.mean", "score.subset.mean", "score.all.samples.mean", "distance.mean", "local_clustering_coefficient.mean", "degree_centrality.mean", "betweenness_centrality.mean", "type_interaction"))))
ranges_df = data.frame(left.range=seq(-1, 0.9, 0.1), right.range=seq(-0.9, 1, 0.1))
for(range.id in rownames(ranges_df)){
  left.range = ranges_df[range.id,"left.range"]
  right.range = ranges_df[range.id,"right.range"]
  subset_df = ppi_gene_coexpression_merged_df %>% filter((score.subset.network >= left.range) & (score.subset.network < right.range))
  subset_direct_df = subset_df %>% filter(distances < 2)
  subset_indirect_df = subset_df %>% filter(distances >= 2)
  scores_df = rbind(scores_df, data.frame(left.range=left.range, right.range=right.range, 
                                          number.of.edges=length(subset_df$difference),
                                          difference.mean=mean(subset_df$difference),
                                          score.subset.mean=mean(subset_df$score.subset.network),
                                          score.all.samples.mean=mean(subset_df$score.all.samples.network),
                                          distance.mean=mean(subset_df$distances),
                                          local_clustering_coefficient.mean=mean(subset_df$local_clustering_coefficient.mean),
                                          degree_centrality.mean=mean(subset_df$degree_centrality.mean),
                                          betweenness_centrality.mean=mean(subset_df$betweenness_centrality.mean),
                                          type_interaction="all"
  ))
  scores_df = rbind(scores_df, data.frame(left.range=left.range, right.range=right.range, 
                                          number.of.edges=length(subset_direct_df$difference),
                                          difference.mean=mean(subset_direct_df$difference),
                                          score.subset.mean=mean(subset_direct_df$score.subset.network),
                                          score.all.samples.mean=mean(subset_direct_df$score.all.samples.network),
                                          distance.mean=mean(subset_direct_df$distances),
                                          local_clustering_coefficient.mean=mean(subset_direct_df$local_clustering_coefficient.mean),
                                          degree_centrality.mean=mean(subset_direct_df$degree_centrality.mean),
                                          betweenness_centrality.mean=mean(subset_direct_df$betweenness_centrality.mean),
                                          type_interaction="direct"
  ))
  scores_df = rbind(scores_df, data.frame(left.range=left.range, right.range=right.range, 
                                          number.of.edges=length(subset_indirect_df$difference),
                                          difference.mean=mean(subset_indirect_df$difference),
                                          score.subset.mean=mean(subset_indirect_df$score.subset.network),
                                          score.all.samples.mean=mean(subset_indirect_df$score.all.samples.network),
                                          distance.mean=mean(subset_indirect_df$distances),
                                          local_clustering_coefficient.mean=mean(subset_indirect_df$local_clustering_coefficient.mean),
                                          degree_centrality.mean=mean(subset_indirect_df$degree_centrality.mean),
                                          betweenness_centrality.mean=mean(subset_indirect_df$betweenness_centrality.mean),
                                          type_interaction="indirect"
  ))
}
rm(subset_df)
rm(subset_direct_df)
rm(subset_indirect_df)


# Create output file
scores_df %>% fwrite(scores_file)


# Get difference of nodes and edges by using different co-expression score thresholds
threshold_results_df = data.frame(matrix(ncol=7,nrow=0, dimnames=list(NULL, c("threshold", "number.of.edges.subset", "number.of.edges.all.samples", "number.of.edges.overlapped", "number.of.edges.lost", "number.of.edges.gained", "type_interaction"))))
thresholds = seq(0.1, 0.9, 0.1)
for (threshold in thresholds){
  network_subset_df = ppi_gene_coexpression_merged_df %>% filter(score.subset.network >= threshold)
  network_subset_df$edge = paste(network_subset_df$Node.1, network_subset_df$Node.2)
  network_all_samples_df = ppi_gene_coexpression_merged_df %>% filter(score.all.samples.network >= threshold)
  network_all_samples_df$edge = paste(network_all_samples_df$Node.1, network_all_samples_df$Node.2)
  edge_overlap = network_subset_df %>% filter(edge %in% network_all_samples_df$edge)
  edges_lost = network_all_samples_df %>% filter(!(edge %in% network_subset_df$edge))
  edges_gained = network_subset_df %>% filter(!(edge %in% network_all_samples_df$edge))
  threshold_results_df = rbind(threshold_results_df, data.frame(threshold=threshold,
                                                                number.of.edges.subset=length(network_subset_df$edge),
                                                                number.of.edges.all.samples=length(network_all_samples_df$edge),
                                                                number.of.edges.overlapped=length(edge_overlap$edge),
                                                                number.of.edges.lost=length(edges_lost$edge),
                                                                number.of.edges.gained=length(edges_gained$edge),
                                                                type_interaction="all"
  ))
  network_subset_direct_df = network_subset_df %>% filter(distances < 2)
  network_all_samples_direct_df = network_all_samples_df %>% filter(distances < 2)
  edge_overlap_direct = edge_overlap %>% filter(distances < 2)
  edges_lost_direct = edges_lost %>% filter(distances < 2)
  edges_gained_direct = edges_gained %>% filter(distances < 2)
  threshold_results_df = rbind(threshold_results_df, data.frame(threshold=threshold, 
                                                                number.of.edges.subset=length(network_subset_direct_df$edge),
                                                                number.of.edges.all.samples=length(network_all_samples_direct_df$edge),
                                                                number.of.edges.overlapped=length(edge_overlap_direct$edge),
                                                                number.of.edges.lost=length(edges_lost_direct$edge),
                                                                number.of.edges.gained=length(edges_gained_direct$edge),
                                                                type_interaction="direct"
  ))
  network_subset_indirect_df = network_subset_df %>% filter(distances >= 2)
  network_all_samples_indirect_df = network_all_samples_df %>% filter(distances >= 2)
  edge_overlap_indirect = edge_overlap %>% filter(distances >= 2)
  edges_lost_indirect = edges_lost %>% filter(distances >= 2)
  edges_gained_indirect = edges_gained %>% filter(distances >= 2)
  threshold_results_df = rbind(threshold_results_df, data.frame(threshold=threshold,
                                                                number.of.edges.subset=length(network_subset_indirect_df$edge),
                                                                number.of.edges.all.samples=length(network_all_samples_indirect_df$edge),
                                                                number.of.edges.overlapped=length(edge_overlap_indirect$edge),
                                                                number.of.edges.lost=length(edges_lost_indirect$edge),
                                                                number.of.edges.gained=length(edges_gained_indirect$edge),
                                                                type_interaction="indirect"
  ))
  
}
rm(network_subset_df)
rm(network_subset_direct_df)
rm(network_subset_indirect_df)
rm(network_all_samples_df)
rm(network_all_samples_direct_df)
rm(network_all_samples_indirect_df)


# Create output file
threshold_results_df %>% fwrite(threshold_results_file)


# Calculate stability depending on the disease genes
disease_results_df = data.frame(matrix(ncol=10,nrow=0, dimnames=list(NULL, c("disease","disease_class","disease_genes_in_network", "num.edges", "distances", "local_clustering_coefficient.mean", "degree_centrality.mean", "betweenness_centrality.mean", "score.subset.network", "type_interaction"))))
for (disease in unique(disease_genes_df$disease)){
  disease_selected_df = disease_genes_df %>% filter(disease == !!disease)
  disease_selected_genes = unique(disease_selected_df$gene_symbol)
  if (length(disease_selected_genes) > 19){
    disease_selected_classes = unique(disease_selected_df$root.concept)
    disease_selected_coexpression_info_df = ppi_gene_coexpression_merged_df %>% filter((Node.1 %in% disease_selected_genes) & (Node.2 %in% disease_selected_genes))
    mean_coexpression_info_all_df = colMeans((disease_selected_coexpression_info_df %>% select(!c(Node.1, Node.2)))) %>% as.data.frame.list()
    coexpression_info_direct_df = disease_selected_coexpression_info_df %>% filter(distances < 2)
    mean_coexpression_info_direct_df = colMeans((coexpression_info_direct_df %>% select(!c(Node.1, Node.2)))) %>% as.data.frame.list()
    coexpression_info_indirect_df = disease_selected_coexpression_info_df %>% filter(distances >= 2)
    mean_coexpression_info_indirect_df = colMeans((coexpression_info_indirect_df %>% select(!c(Node.1, Node.2)))) %>% as.data.frame.list()
    for (disease_selected_class in disease_selected_classes){
      disease_results_df = rbind(disease_results_df, cbind(data.frame(disease=disease, disease_class=disease_selected_class, disease_genes_in_network=length(disease_selected_genes), num.edges=nrow(disease_selected_coexpression_info_df)), mean_coexpression_info_all_df, data.frame(type_interaction="all")))
      disease_results_df = rbind(disease_results_df, cbind(data.frame(disease=disease, disease_class=disease_selected_class, disease_genes_in_network=length(disease_selected_genes), num.edges=nrow(coexpression_info_direct_df)), mean_coexpression_info_direct_df, data.frame(type_interaction="direct")))
      disease_results_df = rbind(disease_results_df, cbind(data.frame(disease=disease, disease_class=disease_selected_class, disease_genes_in_network=length(disease_selected_genes), num.edges=nrow(coexpression_info_indirect_df)), mean_coexpression_info_indirect_df, data.frame(type_interaction="indirect")))
    }
  }
}

# Create output file
disease_results_df %>% fwrite(disease_results_file)

