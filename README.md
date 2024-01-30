# Sample size value of gene expression data in Network Medicine


## Introduction

The scripts folder contains the scripts necessary to run the analysis.

## Code

### 1. Extract the gene expression data

==> folder: extract_data

This is a detailed explanation on how to extract the datasets used in the analysis.

#### GTEx



#### GEO GSE193677

Notebook to run:

```
GEO_data_extraction_GSE193677.Rmd
```

#### TCGA

Check README_TCGA

Notebook to run:

```
TCGA_preparation.Rmd
```

### 2. Process the gene expression datasets

==> folder: process_gene_expression

Notebooks to run:

* For GTEx:

```
GTEx_preprocessing_for_coexpr.Rmd
```

* For GEO GSE193677:

```
GEO_preprocessing_for_coexpr_GSE193677.Rmd
```

* For TCGA:

```
TCGA_preprocessing_for_coexpr.Rmd
```

### 3. Bootstrapping

==> folder: bootstrapping

We create groups of samples of different sizes using random sampling with replacement, starting from size 20 and increasing the size by 20 until reaching the maximum number of samples in the dataset.

For this, we have to run a R markdown notebook for each dataset.

* For GTEx:

```
GTEx_subsampling.Rmd
```

* For GEO GSE193677

```
GEO_subsampling_GSE193677.Rmd
```

* For TCGA:

```
TCGA_subsampling.Rmd
```


### 4. Create the gene co-expression networks

==> folder: gene_coexpression_networks

#### Individual co-expression networks

Execution:

```
Rscript create_gene_coexpression_network_from_samples_list.R -s <samples_file> -f <rnaseq_file> -o <output_file> -m <metric> -n <wto_n> -d <wto_delta> -p 6 -t signed -e pearson -a 0 -c bonferroni
```

This is done for all datasets.

#### Consensus co-expression networks

Then, we calculated the consensus networks between the 5 replicate genes co-expression networks from each sample size and condition.

Execution:

```
Rscript create_consensus_gene_coexpression_network.R -l <list_coexpression_networks_file> -n <consensus_network_file> -t <threshold> -m <method>
```

### 5. Analyze the gene co-expression networks

==> folder: gene_coexpression_networks

Execution:

```
Rscript analyze_coexpression_network_by_significant_edges.R -c <coexpression_network_file> -o <output_results_dir> -s <output_subgraphs_dir> -f <file_name> -t <threshold> -p <ppi_file> -d <disease_genes_file> -e <essential_genes_file> -g <genes_dataset_file>
```


### 6. Parse the results of the analyses and create summary tables

==> folder: parse_results

1. Run: "parse_coexpression_analysis_results.Rmd"

2. Run: "server_example.R"

3. Run: "calculate_rnaseq_datasets_variation.Rmd"

4???? Is necessary?. Run: "calculate_convergence_correlation_types.Rmd"


### 7. Differential co-expression analysis

==> folder: differential_coexpression

#### Calculate differential co-expression analysis

Requires a lot of computational memory: 100000-120000 MB (150000 for a couple of them).

Execution:

```
calculate_differentially_coexpressed_genes.R -d <coexpression_network_file_D> -n <coexpression_network_file_N> -l <output_edges_file> -v <output_nodes_file> -t <threshold> -s <stretch_normalization> -f <filter_by_common_nodes>
```

#### Analyze the results from the differential co-expression analysis

```
Rscript analyze_differentially_coexpressed_genes.R -a <input_dir> -b <name_disease> -c <name_normal> -d <disease_gene_associations_file> -e <disease_name_in_associations_file> -f <ppi_file> -g <ppi_distances_file> -i <plots_dir> -j <tables_dir> -k <pval_threshold> -l <pval_correction> -m <nodes_to_follow_file> -n <drug_targets_file>
```

### 8. Analyze the gene co-expression variation in individual links

==> folder: analyze_coexpression_variation

#### Analysis of gene co-expression variation across five links

We illustrate the variation in correlation values of five randomly selected links with different correlation strengths from GTEx whole blood networks across five replicates of the same sample size (Figure 3B).

Execution:

```
Rscript analysis_specific_edges.R -n <networks_dir> -o <output_file> -m <method> -l <num_links_selected>
```

#### Analysis of standard deviation of links across repetitions of each sample size

First, we merge the consensus networks from different sample sizes into a unique file. This operation requires a lot of computational memory:

```
Rscript merge_networks.R -n <networks_dir> -m <method> -o <output_dataframe> -s <step> -r <max_rep>
```

Then, we read the file and plot the standard deviation of the repetitions across  sample size using a half violin plot (Figure 3C):

```
Rscript analysis_sd_coexpression_weight.R -n <networks_file> -t <output_sd_table_file> -p <output_sd_plot_file>
```

### 9.