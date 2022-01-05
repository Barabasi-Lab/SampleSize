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

# Merge the two dataframes (gene co-expression network and distances network) as both have all possible combinations of pairs of nodes
ppi_gene_coexpression_merged_df <- inner_join(ppi_tissue_sex_filt_distances, tissue_sex_network_merged_df, by=c("Node.1" = "Node.1", "Node.2" = "Node.2"))
rm(ppi_tissue_sex_filt_distances)
rm(tissue_sex_network_merged_df)

# Calculate difference between networks
ppi_gene_coexpression_merged_df$difference = ppi_gene_coexpression_merged_df$score.all.samples.network - ppi_gene_coexpression_merged_df$score.subset.network

difference_df = data.frame(matrix(ncol=2,nrow=0, dimnames=list(NULL, c("difference.range", "number.of.edges"))))
difference_df = rbind(difference_df, data.frame(difference.range="0,0.001", number.of.edges=length((ppi_gene_coexpression_merged_df %>% filter((difference >= 0) & (difference < 0.001)))$difference)))
difference_df = rbind(difference_df, data.frame(difference.range="0.001,0.005", number.of.edges=length((ppi_gene_coexpression_merged_df %>% filter((difference >= 0.001) & (difference < 0.005)))$difference)))
difference_df = rbind(difference_df, data.frame(difference.range="0.005,0.01", number.of.edges=length((ppi_gene_coexpression_merged_df %>% filter((difference >= 0.005) & (difference < 0.01)))$difference)))
difference_df = rbind(difference_df, data.frame(difference.range="0.01,0.05", number.of.edges=length((ppi_gene_coexpression_merged_df %>% filter((difference >= 0.01) & (difference < 0.05)))$difference)))
difference_df = rbind(difference_df, data.frame(difference.range="0.05,0.1", number.of.edges=length((ppi_gene_coexpression_merged_df %>% filter((difference >= 0.05) & (difference < 0.1)))$difference)))
difference_df = rbind(difference_df, data.frame(difference.range="0.1", number.of.edges=length((ppi_gene_coexpression_merged_df %>% filter(difference >= 0.1))$difference)))
difference_df = rbind(difference_df, data.frame(difference.range="-0.001,0", number.of.edges=length((ppi_gene_coexpression_merged_df %>% filter((difference < 0) & (difference >= -0.001)))$difference)))
difference_df = rbind(difference_df, data.frame(difference.range="-0.005,-0.001", number.of.edges=length((ppi_gene_coexpression_merged_df %>% filter((difference < -0.001) & (difference >= -0.005)))$difference)))
difference_df = rbind(difference_df, data.frame(difference.range="-0.01,-0.005", number.of.edges=length((ppi_gene_coexpression_merged_df %>% filter((difference < -0.005) & (difference >= -0.01)))$difference)))
difference_df = rbind(difference_df, data.frame(difference.range="-0.05,-0.01", number.of.edges=length((ppi_gene_coexpression_merged_df %>% filter((difference < -0.01) & (difference >= -0.05)))$difference)))
difference_df = rbind(difference_df, data.frame(difference.range="-0.1,-0.05", number.of.edges=length((ppi_gene_coexpression_merged_df %>% filter((difference < -0.05) & (difference >= -0.1)))$difference)))
difference_df = rbind(difference_df, data.frame(difference.range="-0.1", number.of.edges=length((ppi_gene_coexpression_merged_df %>% filter(difference < -0.1))$difference)))

# Define positives and negatives
positives = ppi_gene_coexpression_merged_df %>% filter(distances < 2)
negatives = ppi_gene_coexpression_merged_df %>% filter(distances > 4)

# Calculate a discrete score based on the quantiles
ppi_gene_coexpression_merged_df$score_discrete = ifelse( ppi_gene_coexpression_merged_df$score >= quantiles_score[[4]], '75-100%', ifelse( (ppi_gene_coexpression_merged_df$score >= quantiles_score[[3]]) & (ppi_gene_coexpression_merged_df$score < quantiles_score[[4]]), '50-75%', ifelse( (ppi_gene_coexpression_merged_df$score >= quantiles_score[[2]]) & (ppi_gene_coexpression_merged_df$score < quantiles_score[[3]]), '25-50%', '0-25%') ) )
# Create a table comparing score to distance
table_score_distance = table(ppi_gene_coexpression_merged_df$score_discrete, ppi_gene_coexpression_merged_df$distances)

# Calculate general performance on predicting distance based on the whole dataset
general_pred_df = data.frame(TP=(nrow(positives %>% filter(score >= quantiles_score[[3]]))),
                             FP=(nrow(positives %>% filter(score < quantiles_score[[3]]))),
                             TN=(nrow(negatives %>% filter(score >= quantiles_score[[3]]))),
                             FN=(nrow(negatives %>% filter(score < quantiles_score[[3]])))
)
general_pred_df$TPR = ifelse((general_pred_df$TP + general_pred_df$FN) > 0, general_pred_df$TP / (general_pred_df$TP + general_pred_df$FN), 0)
general_pred_df$FPR = general_pred_df$FP / (general_pred_df$FP + general_pred_df$TN)
general_pred_df$accuracy = (general_pred_df$TP + general_pred_df$TN) / (general_pred_df$TP + general_pred_df$TN + general_pred_df$FP + general_pred_df$FN)
general_pred_df$F1 = 2*general_pred_df$TP / (2*general_pred_df$TP + general_pred_df$FP + general_pred_df$FN)
general_pred_df$MCC = mcc(TP=general_pred_df$TP, FP=general_pred_df$FP, TN=general_pred_df$TN, FN=general_pred_df$FN)
general_pred_df$corr = corr_distance_score
general_pred_df$type_prediction = 'general'

# Calculate prediction of distances with expression based on balanced subsets
pred_df = data.frame(matrix(ncol=11,nrow=0, dimnames=list(NULL, c("rep", "TP", "FP", "TN", "FN", "TPR", "FPR", "accuracy", "F1", "MCC", "corr"))))
for (r in 1:1000){
  #print(r)
  positive_set = sample(rownames(positives), size=ifelse(length(rownames(positives)) < 100, length(rownames(positives)), 100), replace=FALSE)
  negative_set = sample(rownames(negatives), size=ifelse(length(rownames(negatives)) < 100, length(rownames(negatives)), 100), replace=FALSE)
  positive_set = positives %>% slice(as.numeric(positive_set))
  negative_set = negatives %>% slice(as.numeric(negative_set))
  testing_set = rbind(positive_set, negative_set)
  corr_testing_set = cor(x=testing_set$distances, y=testing_set$score, method="pearson")
  TP = nrow(testing_set %>% filter((distances < 2) & (score >= quantiles_score[[3]])))
  FP = nrow(testing_set %>% filter((distances < 2) & (score < quantiles_score[[3]])))
  TN = nrow(testing_set %>% filter((distances > 4) & (score < quantiles_score[[3]])))
  FN = nrow(testing_set %>% filter((distances > 4) & (score >= quantiles_score[[3]])))
  TPR = ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  FPR = FP / (FP + TN)
  accuracy = (TP + TN) / (TP + TN + FP + FN)
  F1 = 2*TP / (2*TP + FP + FN)
  MCC = mcc(TP=TP, FP=FP, TN=TN, FN=FN)
  pred_df = rbind(pred_df, data.frame(rep=r, TP=TP, FP=FP, TN=TN, FN=FN, TPR=TPR, FPR=FPR, accuracy=accuracy, F1=F1, MCC=MCC, corr=corr_testing_set))
}
mean_pred = pred_df %>% colMeans() %>% as.data.frame.list() %>% select(!(rep))
mean_pred$type_prediction = 'balanced'
pred = rbind(general_pred_df, mean_pred)


# Calculate prediction of high clustering PPIs based on co-expression
positives_clustering = ppi_gene_coexpression_merged_df %>% filter((Node.1 %in% V(high_clustering_graph)$name) & (Node.2 %in% V(high_clustering_graph)$name) & (distances < 2))
negatives_clustering_direct = ppi_gene_coexpression_merged_df %>% filter((Node.1 %in% V(low_clustering_graph)$name) & (Node.2 %in% V(low_clustering_graph)$name) & (distances < 2))
negatives_clustering_indirect = ppi_gene_coexpression_merged_df %>% filter(distances > 4)
negatives_clustering = negatives_clustering_indirect
positive_clustering_score_dist = data.frame(Node.1=positives_clustering$Node.1, Node.2=positives_clustering$Node.2, score=positives_clustering$score, type="positive")
negative_clustering_score_dist = data.frame(Node.1=negatives_clustering$Node.1, Node.2=negatives_clustering$Node.2, score=negatives_clustering$score, type="negative")
clustering_score_dist = rbind(positive_clustering_score_dist, negative_clustering_score_dist)
#clustering_score_dist %>%
#  ggplot( aes(x=score, fill=type)) +
#  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
#  scale_fill_manual(values=c("#69b3a2", "#404080"))

positive_clustering_nodes = unique(c(as.vector(positives_clustering$Node.1), as.vector(positives_clustering$Node.2)))
negative_clustering_nodes = unique(c(as.vector(negatives_clustering$Node.1), as.vector(negatives_clustering$Node.2)))
clustering_pred_df = data.frame(matrix(ncol=11,nrow=0, dimnames=list(NULL, c("rep", "TP", "FP", "TN", "FN", "TPR", "FPR", "accuracy", "F1", "MCC", "corr"))))
for (r in 1:1000){
  #print(r)
  positive_set = sample(rownames(positives_clustering), size=ifelse(length(rownames(positives_clustering)) < 100, length(rownames(positives_clustering)), 100), replace=FALSE)
  negative_set = sample(rownames(negatives_clustering), size=ifelse(length(rownames(negatives_clustering)) < 100, length(rownames(negatives_clustering)), 100), replace=FALSE)
  positive_set = positives_clustering %>% slice(as.numeric(positive_set))
  negative_set = negatives_clustering %>% slice(as.numeric(negative_set))
  testing_set = rbind(positive_set, negative_set)
  corr_testing_set = cor(x=testing_set$distances, y=testing_set$score, method="pearson")
  TP = nrow(testing_set %>% filter((Node.1 %in% positive_clustering_nodes) & (Node.2 %in% positive_clustering_nodes) & (score >= quantiles_score[[3]])))
  FP = nrow(testing_set %>% filter((Node.1 %in% positive_clustering_nodes) & (Node.2 %in% positive_clustering_nodes) & (score < quantiles_score[[3]])))
  TN = nrow(testing_set %>% filter((Node.1 %in% negative_clustering_nodes) & (Node.2 %in% negative_clustering_nodes) & (score < quantiles_score[[3]])))
  FN = nrow(testing_set %>% filter((Node.1 %in% negative_clustering_nodes) & (Node.2 %in% negative_clustering_nodes) & (score >= quantiles_score[[3]])))
  TPR = ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  FPR = FP / (FP + TN)
  accuracy = (TP + TN) / (TP + TN + FP + FN)
  F1 = 2*TP / (2*TP + FP + FN)
  MCC = mcc(TP=TP, FP=FP, TN=TN, FN=FN)
  clustering_pred_df = rbind(clustering_pred_df, data.frame(rep=r, TP=TP, FP=FP, TN=TN, FN=FN, TPR=TPR, FPR=FPR, accuracy=accuracy, F1=F1, MCC=MCC, corr=corr_testing_set))
}
clustering_mean_pred = clustering_pred_df %>% colMeans() %>% as.data.frame.list() %>% select(!(rep))
clustering_mean_pred$type_prediction = 'clustering'
pred = rbind(pred, clustering_mean_pred)


# Calculate prediction of high degree PPIs based on co-expression
positives_degree = ppi_gene_coexpression_merged_df %>% filter((Node.1 %in% V(high_degree_graph)$name) & (Node.2 %in% V(high_degree_graph)$name) & (distances < 2))
negatives_degree_direct = ppi_gene_coexpression_merged_df %>% filter((Node.1 %in% V(low_degree_graph)$name) & (Node.2 %in% V(low_degree_graph)$name) & (distances < 2))
negatives_degree_indirect = ppi_gene_coexpression_merged_df %>% filter(distances > 4)
negatives_degree = negatives_degree_indirect
positive_degree_score_dist = data.frame(Node.1=positives_degree$Node.1, Node.2=positives_degree$Node.2, score=positives_degree$score, type="positive")
negative_degree_score_dist = data.frame(Node.1=negatives_degree$Node.1, Node.2=negatives_degree$Node.2, score=negatives_degree$score, type="negative")
degree_score_dist = rbind(positive_degree_score_dist, negative_degree_score_dist)
#degree_score_dist %>%
#  ggplot( aes(x=score, fill=type)) +
#  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
#  scale_fill_manual(values=c("#69b3a2", "#404080"))

positive_degree_nodes = unique(c(as.vector(positives_degree$Node.1), as.vector(positives_degree$Node.2)))
negative_degree_nodes = unique(c(as.vector(negatives_degree$Node.1), as.vector(negatives_degree$Node.2)))
degree_pred_df = data.frame(matrix(ncol=11,nrow=0, dimnames=list(NULL, c("rep", "TP", "FP", "TN", "FN", "TPR", "FPR", "accuracy", "F1", "MCC", "corr"))))
for (r in 1:1000){
  #print(r)
  positive_set = sample(rownames(positives_degree), size=ifelse(length(rownames(positives_degree)) < 100, length(rownames(positives_degree)), 100), replace=FALSE)
  negative_set = sample(rownames(negatives_degree), size=ifelse(length(rownames(negatives_degree)) < 100, length(rownames(negatives_degree)), 100), replace=FALSE)
  positive_set = positives_degree %>% slice(as.numeric(positive_set))
  negative_set = negatives_degree %>% slice(as.numeric(negative_set))
  testing_set = rbind(positive_set, negative_set)
  corr_testing_set = cor(x=testing_set$distances, y=testing_set$score, method="pearson")
  TP = nrow(testing_set %>% filter((Node.1 %in% positive_degree_nodes) & (Node.2 %in% positive_degree_nodes) & (score >= quantiles_score[[3]])))
  FP = nrow(testing_set %>% filter((Node.1 %in% positive_degree_nodes) & (Node.2 %in% positive_degree_nodes) & (score < quantiles_score[[3]])))
  TN = nrow(testing_set %>% filter((Node.1 %in% negative_degree_nodes) & (Node.2 %in% negative_degree_nodes) & (score < quantiles_score[[3]])))
  FN = nrow(testing_set %>% filter((Node.1 %in% negative_degree_nodes) & (Node.2 %in% negative_degree_nodes) & (score >= quantiles_score[[3]])))
  TPR = ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  FPR = FP / (FP + TN)
  accuracy = (TP + TN) / (TP + TN + FP + FN)
  F1 = 2*TP / (2*TP + FP + FN)
  MCC = mcc(TP=TP, FP=FP, TN=TN, FN=FN)
  degree_pred_df = rbind(degree_pred_df, data.frame(rep=r, TP=TP, FP=FP, TN=TN, FN=FN, TPR=TPR, FPR=FPR, accuracy=accuracy, F1=F1, MCC=MCC, corr=corr_testing_set))
}
degree_mean_pred = degree_pred_df %>% colMeans() %>% as.data.frame.list() %>% select(!(rep))
degree_mean_pred$type_prediction = 'degree'
pred = rbind(pred, degree_mean_pred)


# Calculate prediction of high betweenness PPIs based on co-expression
positive_betweenness = ppi_gene_coexpression_merged_df %>% filter((Node.1 %in% V(high_betweenness_graph)$name) & (Node.2 %in% V(high_betweenness_graph)$name) & (distances < 2))
#negative_betweenness_direct = ppi_gene_coexpression_merged_df %>% filter((Node.1 %in% V(low_betweenness_graph)$name) & (Node.2 %in% V(low_betweenness_graph)$name) & (distances < 2))
negative_betweenness_indirect = ppi_gene_coexpression_merged_df %>% filter(distances > 4)
negative_betweenness = negative_betweenness_indirect
positive_betweenness_score_dist = data.frame(Node.1=positive_betweenness$Node.1, Node.2=positive_betweenness$Node.2, score=positive_betweenness$score, type="positive")
negative_betweenness_score_dist = data.frame(Node.1=negative_betweenness$Node.1, Node.2=negative_betweenness$Node.2, score=negative_betweenness$score, type="negative")
betweenness_score_dist = rbind(positive_betweenness_score_dist, negative_betweenness_score_dist)
#betweenness_score_dist %>%
#  ggplot( aes(x=score, fill=type)) +
#  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
#  scale_fill_manual(values=c("#69b3a2", "#404080"))

positive_betweenness_nodes = unique(c(as.vector(positive_betweenness$Node.1), as.vector(positive_betweenness$Node.2)))
negative_betweenness_nodes = unique(c(as.vector(negative_betweenness$Node.1), as.vector(negative_betweenness$Node.2)))
betweenness_pred_df = data.frame(matrix(ncol=11,nrow=0, dimnames=list(NULL, c("rep", "TP", "FP", "TN", "FN", "TPR", "FPR", "accuracy", "F1", "MCC", "corr"))))
for (r in 1:1000){
  #print(r)
  positive_set = sample(rownames(positive_betweenness), size=ifelse(length(rownames(positive_betweenness)) < 100, length(rownames(positive_betweenness)), 100), replace=FALSE)
  negative_set = sample(rownames(negative_betweenness), size=ifelse(length(rownames(negative_betweenness)) < 100, length(rownames(negative_betweenness)), 100), replace=FALSE)
  positive_set = positive_betweenness %>% slice(as.numeric(positive_set))
  negative_set = negative_betweenness %>% slice(as.numeric(negative_set))
  testing_set = rbind(positive_set, negative_set)
  corr_testing_set = cor(x=testing_set$distances, y=testing_set$score, method="pearson")
  TP = nrow(testing_set %>% filter((Node.1 %in% positive_betweenness_nodes) & (Node.2 %in% positive_betweenness_nodes) & (score >= quantiles_score[[3]])))
  FP = nrow(testing_set %>% filter((Node.1 %in% positive_betweenness_nodes) & (Node.2 %in% positive_betweenness_nodes) & (score < quantiles_score[[3]])))
  TN = nrow(testing_set %>% filter((Node.1 %in% negative_betweenness_nodes) & (Node.2 %in% negative_betweenness_nodes) & (score < quantiles_score[[3]])))
  FN = nrow(testing_set %>% filter((Node.1 %in% negative_betweenness_nodes) & (Node.2 %in% negative_betweenness_nodes) & (score >= quantiles_score[[3]])))
  TPR = ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  FPR = FP / (FP + TN)
  accuracy = (TP + TN) / (TP + TN + FP + FN)
  F1 = 2*TP / (2*TP + FP + FN)
  MCC = mcc(TP=TP, FP=FP, TN=TN, FN=FN)
  betweenness_pred_df = rbind(betweenness_pred_df, data.frame(rep=r, TP=TP, FP=FP, TN=TN, FN=FN, TPR=TPR, FPR=FPR, accuracy=accuracy, F1=F1, MCC=MCC, corr=corr_testing_set))
}
betweenness_mean_pred = betweenness_pred_df %>% colMeans() %>% as.data.frame.list() %>% select(!(rep))
betweenness_mean_pred$type_prediction = 'betweenness'
pred = rbind(pred, betweenness_mean_pred)

# Create output file
#output_file = paste(output_dir, "comparison_coexpression_ppi.txt", sep='/')
pred %>% fwrite(output_file)




