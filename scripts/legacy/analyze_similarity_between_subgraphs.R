#!/usr/bin/env Rscript
packrat::init("/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeR")
library(optparse)
library(data.table)
library(dplyr)
library(ggplot2)
library(igraph)
library(purrr)
set.seed(1510)



graph_similarity = function(g1, g2){
  # compute the Jaccard and overlap Index
  
  
  g_intersect = graph.intersection(g1, g2, keep.all.vertices = FALSE)
  g_union = graph.union(g1, g2, byname = T)
  
  #Jaccard index
  jAB = ecount(g_intersect)/ecount(g_union)
  
  # Overlap index
  oAB = ecount(g_intersect)/min(ecount(g1), ecount(g2))
  return(c(jaccardIndex = jAB, overlapindex = oAB))
}



selected_diseases = c("arthritis, rheumatoid")
selected_diseases <- gsub(' ', '.', gsub(', ', '.', gsub('-', '.', gsub('[\\(\\)]', '', selected_diseases))))


# Get network files
networks_info_df = data.frame(matrix(ncol=10,nrow=0, dimnames=list(NULL, c("file", "method", "dataset", "type_dataset", "type_analysis", "tissue", "sex", "size", "rep", "edge_threshold"))))
input_dir = "/home/j.aguirreplans/Projects/Scipher/SampleSize/data/out"
for(dataset in c('scipher', 'gtex')){
  types_dataset = list.files(paste(input_dir, paste("networks_", dataset, sep=""), sep="/"))
  for(type_dataset in types_dataset){
    if (dataset == "gtex"){
      type_dataset_split = strsplit(type_dataset, split="_")[[1]]
      if (length(type_dataset_split == 2)){
        tissue = type_dataset_split[1]
        sex = type_dataset_split[2]
      } else {
        tissue = type_dataset_split[1]
        sex = "both"
      }
    } else {
      tissue = NA
      sex = NA
    }
    for(type_analysis in c("subgraphs", "main_core", "essential_genes", "disease_genes")){
      network_files = list.files(paste(input_dir, paste("networks_", dataset, sep=""), type_dataset, type_analysis, sep="/"))
      for (network_file in network_files){
        network_file_path = paste(input_dir, paste("networks_", dataset, sep=""), type_dataset, type_analysis, network_file, sep="/")
        network_file_proc = gsub("_main_core_subgraph.txt|_essential_genes_subgraph.txt|_subgraph.txt|_disease_genes_subgraphs","",network_file)
        file_split = strsplit(network_file_proc, split="_")[[1]]
        method = file_split[1]
        edge_threshold = as.numeric(file_split[(length(file_split))])
        size = file_split[(length(file_split)-4)]
        r = file_split[(length(file_split))-2]
        if(r=="samples"){
          size = "all"
          r = NA
        } else {
          size = as.numeric(size)
          r = as.numeric(r)
        }
        if(type_analysis == "disease_genes"){
          networks_info_df <- rbind(networks_info_df, data.frame(file=paste(network_file_path, paste(selected_diseases, ".txt", sep=""), sep="/"), method=method, dataset=dataset, type_dataset=type_dataset, type_analysis=type_analysis, tissue=tissue, sex=sex, size=size, rep=r, edge_threshold=edge_threshold))
        } else {
          networks_info_df <- rbind(networks_info_df, data.frame(file=network_file_path, method=method, dataset=dataset, type_dataset=type_dataset, type_analysis=type_analysis, tissue=tissue, sex=sex, size=size, rep=r, edge_threshold=edge_threshold))
        }
      }
    }
  }
}



# Analyze groups of network files by size
output_dir = '/home/j.aguirreplans/Projects/Scipher/SampleSize/data/data_shiny_app'
output_network_similarity_file = paste(output_dir, 'network_similarity_by_size.csv', sep='/')
network_similarity_df = data.frame(matrix(ncol=11,nrow=0, dimnames=list(NULL, c("file1", "file2", "dataset", "type_dataset", "type_analysis", "tissue", "sex", "size", "edge_threshold", "jaccardIndex", "overlapindex"))))
edge_threshold = 100000

if(!(file.exists(output_network_similarity_file))){
  for(method in c("spearman", "wgcna")){
    print(method)
    for(dataset in c("scipher", "gtex")){
      print(dataset)
      types_dataset = list.files(paste(input_dir, paste("networks_", dataset, sep=""), sep="/"))
      for(type_dataset in types_dataset){
        if (dataset == "gtex"){
          type_dataset_split = strsplit(type_dataset, split="_")[[1]]
          if (length(type_dataset_split == 2)){
            tissue = type_dataset_split[1]
            sex = type_dataset_split[2]
          } else {
            tissue = type_dataset_split[1]
            sex = "both"
          }
        } else {
          tissue = NA
          sex = NA
        }
        for(type_analysis in c("subgraphs", "main_core", "essential_genes", "disease_genes")){
          print(type_analysis)
          selected_networks_df = networks_info_df %>% filter((method == !!method) & (dataset == !!dataset) & (type_dataset == !!type_dataset) & (type_analysis == !!type_analysis) & (edge_threshold == !!edge_threshold))
          sizes = sort(as.numeric(unique(selected_networks_df$size[!(selected_networks_df$size == "all")])), decreasing = T)
          for (size in sizes){
            selected_network_files = sort((selected_networks_df %>% filter(size == !!size))$file)
            pairs_df = t(combn(selected_network_files,2))
            for(sel_row in 1:nrow(pairs_df)){
              pair = pairs_df[sel_row,]
              if((file.exists(pair[1])) & (file.exists(pair[2]))){
                if(type_analysis == "subgraphs"){
                  network1 = graph_from_data_frame(fread(pair[1]))
                  network2 = graph_from_data_frame(fread(pair[2]))
                } else {
                  network1 = read_graph(file = pair[1], format = "graphml")
                  network2 = read_graph(file = pair[2], format = "graphml")
                }
                graph_sim_res = graph_similarity(g1=network1, g2=network2)
                network_similarity_df <- rbind(network_similarity_df, data.frame(file1=pair[1], file2=pair[2], method=method, dataset=dataset, type_dataset=type_dataset, type_analysis=type_analysis, tissue=tissue, sex=sex, size=size, edge_threshold=edge_threshold, jaccardIndex=graph_sim_res[[1]], overlapindex=graph_sim_res[[2]]))
              }
            }
          }
        }
      }
    }
  }
  network_similarity_df %>% fwrite(output_network_similarity_file)s}


#network_similarity_df$size = as.numeric(network_similarity_df$size)
#network_similarity_by_size = network_similarity_df %>%
#  group_by(size) %>%
#  summarise_at(vars(overlapindex), list(mean=mean, median=median, sd=sd, var=var)) %>%
#  arrange(size)
#network_similarity_by_size$mean.upper = network_similarity_by_size$mean + network_similarity_by_size$var
#network_similarity_by_size$mean.lower = network_similarity_by_size$mean - network_similarity_by_size$var

# Figure based on https://www.r-graph-gallery.com/45-confidence-interval-around-polynomial-curve-fitting.html
#plot(network_similarity_df$size, network_similarity_df$overlapindex,col=rgb(0.4,0.4,0.8,0.6),pch=16 , cex=1.3) 
#lines(network_similarity_by_size$size, network_similarity_by_size$mean, col=2, lwd=2 )  
#polygon(c(network_similarity_by_size$size, rev(network_similarity_by_size$size)), c(network_similarity_by_size$mean.upper, rev(network_similarity_by_size$mean.lower)), col = rgb(0.7,0.7,0.7,0.4) , border = NA)



# Compare network from all samples with networks from subsamples
output_network_comparison_all_samples_file = paste(output_dir, 'network_comparison_subsamples_vs_all_samples.csv', sep='/')
network_comparison_all_samples_df = data.frame(matrix(ncol=11,nrow=0, dimnames=list(NULL, c("file", "dataset", "type_dataset", "type_analysis", "tissue", "sex", "size", "rep", "edge_threshold", "jaccardIndex", "overlapindex"))))

if(!(file.exists(output_network_comparison_all_samples_file))){
  for(method in c("spearman", "wgcna")){
    print(method)
    for(dataset in c("scipher", "gtex")){
      print(dataset)
      types_dataset = list.files(paste(input_dir, paste("networks_", dataset, sep=""), sep="/"))
      for(type_dataset in types_dataset){
        if (dataset == "gtex"){
          type_dataset_split = strsplit(type_dataset, split="_")[[1]]
          if (length(type_dataset_split == 2)){
            tissue = type_dataset_split[1]
            sex = type_dataset_split[2]
          } else {
            tissue = type_dataset_split[1]
            sex = "both"
          }
        } else {
          tissue = NA
          sex = NA
        }
        for(type_analysis in c("subgraphs", "main_core", "essential_genes", "disease_genes")){
          print(type_analysis)
          selected_networks_df = networks_info_df %>% filter((method == !!method) & (dataset == !!dataset) & (type_dataset == !!type_dataset) & (type_analysis == !!type_analysis) & (edge_threshold == !!edge_threshold) & (!(size == "all")))
          selected_network_all_samples_df = networks_info_df %>% filter((method == !!method) & (dataset == !!dataset) & (type_dataset == !!type_dataset) & (type_analysis == !!type_analysis) & (edge_threshold == !!edge_threshold) & (size == "all"))
          if(type_analysis == "subgraphs"){
            network_all_samples = graph_from_data_frame(fread(selected_network_all_samples_df$file))
          } else {
            network_all_samples = read_graph(file = selected_network_all_samples_df$file, format = "graphml")
          }
          for(sel_row in 1:nrow(selected_networks_df)){
            sel_net_info = selected_networks_df[sel_row,]
            if(file.exists(sel_net_info$file)){
              if(type_analysis == "subgraphs"){
                network_subsamples = graph_from_data_frame(fread(sel_net_info$file))
              } else {
                network_subsamples = read_graph(file = sel_net_info$file, format = "graphml")
              }
              graph_sim_res = graph_similarity(g1=network_subsamples, g2=network_all_samples)
              network_comparison_all_samples_df <- rbind(network_comparison_all_samples_df, cbind(sel_net_info, data.frame(jaccardIndex=graph_sim_res[[1]], overlapindex=graph_sim_res[[2]])))
            }
          }
        }
      }
    }
  }
  network_comparison_all_samples_df %>% fwrite(output_network_comparison_all_samples_file)
}



#network_comparison_all_samples_df$size = as.numeric(network_comparison_all_samples_df$size)
#network_comparison_all_samples_by_size = network_comparison_all_samples_df %>%
#  group_by(size) %>%
#  summarise_at(vars(overlapindex), list(mean=mean, median=median, sd=sd, var=var)) %>%
#  arrange(size)
#network_comparison_all_samples_by_size$mean.upper = network_comparison_all_samples_by_size$mean + network_comparison_all_samples_by_size$var
#network_comparison_all_samples_by_size$mean.lower = network_comparison_all_samples_by_size$mean - network_comparison_all_samples_by_size$var

#plot(network_comparison_all_samples_df$size, network_comparison_all_samples_df$overlapindex,col=rgb(0.4,0.4,0.8,0.6),pch=16 , cex=1.3) 
#lines(network_comparison_all_samples_by_size$size, network_comparison_all_samples_by_size$mean, col=2, lwd=2 )  
#polygon(c(network_comparison_all_samples_by_size$size, rev(network_comparison_all_samples_by_size$size)), c(network_comparison_all_samples_by_size$mean.upper, rev(network_comparison_all_samples_by_size$mean.lower)), col = rgb(0.7,0.7,0.7,0.4) , border = NA)


