# Enhancing the accuracy of Network Medicine through understanding the impact of sample size in gene co-expression networks

## Description

This repository contains the scripts used to analyze the impact of sample size on gene co-expression networks. We used the following datasets:

* GTEx
* TCGA
* Rheumatoid Arthritis dataset

We created the gene co-expression networks using **Pearson correlation** and employed **bootstrapping** to generate gene co-expression networks of different sample sizes. We also calculated the **consensus networks** between the 5 replicate genes co-expression networks from each sample size and condition.

We made the following gene co-expression network analyses:

* **Link discovery rate analysis**: We analyzed the link discovery rate of gene expression datasets across gene co-expression networks of different sample sizes.
* **Power-law model fitting**: We fitted the gene co-expression networks of each dataset to a power-law model that describes the rate of discovery of statistically significant links across sample size.
* **Co-expression variation analysis**: We analyzed the variation in co-expression values across different network repetitions and sample sizes.
* **Differential co-expression analysis**: We calculated the differential co-expression between disease and normal gene co-expression networks.
* **Protein-protein interaction co-expression analysis**: We analyzed the influence of sample size on protein-protein interaction co-expression.

We also created summary tables and figures to illustrate the results.


## Table of contents

1. [Code](#code)
    1. [Extract the gene expression data](#1-extract-the-gene-expression-data)
    2. [Process the gene expression datasets](#2-process-the-gene-expression-datasets)
    3. [Bootstrapping](#3-bootstrapping)
    4. [Create the gene co-expression networks](#4-create-the-gene-co-expression-networks)
    5. [Analyze the gene co-expression networks](#5-analyze-the-gene-co-expression-networks)
    6. [Parse the results of the analyses and create summary tables](#6-parse-the-results)
    7. [Differential co-expression analysis](#7-differential-co-expression-analysis)
    8. [Analyze the gene co-expression variation in individual links](#8-analyze-the-gene-co-expression-variation-in-individual-links)
    9. [Analyze the gene co-expression of protein-protein interactions](#9-analyze-the-gene-co-expression-of-protein-protein-interactions)
    10. [Create the figures](#10-create-the-figures)


## Code

### 1. Extract the gene expression data

=> folder: `scripts/extract_data`

This is a detailed explanation on how to extract the datasets used in the analysis.

#### 1.1. GTEx

The data is downloaded from the GTEx portal (https://gtexportal.org/home/).

We downloaded the datasets from v8 release:

* **Gene expression reads**: GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct
* **Sample information**: GTEx/v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
* **Subject information**: GTEx/v8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt

We executed the following notebook to analyze the data and prepare it for the analysis:

```
scripts/extract_data/GTEx_analysis_for_coexpr.Rmd
```

#### 1.3. TCGA

##### 1.3.1. Download the GDC-Client App

We downloaded the GDC-CLIENT app (Ubuntu - Client), necessary to download the dataset, at the following website: https://gdc.cancer.gov/access-data/gdc-data-transfer-tool

This video was also useful to get to know how to download the data: https://www.youtube.com/watch?v=GDxj8DrkZok 

##### 1.3.2. Use the GenomicDataCommons R package to download the metadata and manifest

To download the TCGA dataset, we need to define the the parameters that we want to download and create a manifest file that will be used by the GDC-CLIENT app to download the dataset.

To analyze the data from TCGA, decide which parameters we wanted to use and create the manifest file, we used the GenomicDataCommons R package.

The whole process is detailed in the following Rmarkdown script: 

```
create_TCGA_manifest.Rmd
```

The script will generate the following files:
* Manifest file: `gdc_manifest.2022-11-18.txt`
* Metadata file: `metadata.txt`

##### 1.3.3. Download the dataset

To download the dataset, we executed the following command:

```
gdc-client download -m /path/to/Databases/TCGA/2022-11-18-Dataset/TCGA/raw/additional/gdc_manifest.2022-11-18.txt -d /path/to/Databases/TCGA/2022-11-18-Dataset/TCGA/raw/data
```

The process of downloading the dataset requires a lot of computational memory, so we used a computational cluster to run the command.

##### 1.3.4. Compile the dataset

The TCGA downloaded dataset is scattered in many files. To compile it into a unique file, we executed the following command:

```
python scripts/extract_data/compile_tcga.py -d /path/to/Databases/TCGA/2022-11-18-Dataset/TCGA/raw/data -o /path/to/Databases/TCGA/2022-11-18-Dataset/TCGA/out/reads/TCGA_reads.csv -t unstranded
```

The process of compiling the dataset requires a lot of computational memory, so we used a computational cluster to run the command.

##### 1.3.5. Analyze and prepare the dataset

We executed the following notebook to analyze the data and prepare it for the analysis:

```
scripts/extract_data/TCGA_preparation.Rmd
```


### 2. Process the gene expression datasets

=> folder: `scripts/process_gene_expression`

The gene expression datasets need to be pre-processed to prepare them for the generation of gene co-expression networks. We have to transform the gene identifiers to HGNC names, filter out the genes with low expression counts, and format the data into a matrix with genes in rows and samples in columns. To do this, we have to run the following Rmarkdown notebooks:

* For GTEx:

```
scripts/process_gene_expression/GTEx_preprocessing_for_coexpr.Rmd
```

* For the Rheumatoid Arthritis dataset:

```
scripts/process_gene_expression/scipher_preprocessing.Rmd
```

* For TCGA:

```
scripts/process_gene_expression/TCGA_preprocessing_for_coexpr.Rmd
```

### 3. Bootstrapping

=> folder: `scripts/bootstrapping`

We create groups of samples of different sizes using random sampling with replacement, starting from size 20 and increasing the size by 20 until reaching the maximum number of samples in the dataset. For this, we have to run a R markdown notebook for each dataset:

* For GTEx:

```
scripts/bootstrapping/GTEx_subsampling.Rmd
```

* For Rheumatoid Arthritis dataset:

```
scripts/bootstrapping/scipher_subsampling.Rmd
```

* For TCGA:

```
scripts/bootstrapping/TCGA_subsampling.Rmd
```


### 4. Create the gene co-expression networks

=> folder: `scripts/gene_coexpression_networks`

Here, we describe the process to create the gene co-expression networks using the bootstrapped datasets.

#### Individual gene co-expression networks

We calculate the gene co-expression networks using the bootstrapped datasets. To do it, we ran the following command for each subsampled dataset:

```
Rscript scripts/gene_coexpression_networks/create_gene_coexpression_network_from_samples_list.R -s <samples_file> -f <rnaseq_file> -o <output_file> -m <metric> -n <wto_n> -d <wto_delta> -p <wgcna_power> -t <wgcna_type> -e <mi_estimator> -a <aracne_eps> -c <correction_method>
```

Where:
* `samples_file`: file with the list of samples to be used in the analysis.
* `rnaseq_file`: file with the gene expression data.
* `output_file`: file to save the gene co-expression network.
* `metric`: metric to calculate the co-expression of pairs of genes (e.g., pearson, spearman, mutual_information, wto, wgcna, aracne).
* `wto_n`: Number of wTO bootstrap repetitions to calculate the p-value in wto method.
* `wto_delta`: Value that defines the interval of confidence from which the p-values of the bootstrap repetitions are calculated in wto method.
* `wgcna_power`: Power parameter to calculate the adjacency matrix in wgcna method.
* `wgcna_type`: Type of adjacency matrix to calculate in wgcna method.
* `mi_estimator`: Estimator to calculate the mutual information in mutual_information method.
* `aracne_eps`: Epsilon parameter to calculate the mutual information in aracne method.
* `correction_method`: Method to correct the p-values of the co-expression values (e.g., bonferroni, fdr).

By default, we use as metrix `pearson` and as correction_method `bonferroni`.
Here we show an example of execution with the GTEx dataset for sample size 100 and repetition 1:

```
Rscript scripts/gene_coexpression_networks/create_gene_coexpression_network_from_samples_list.R -s /path/to/GTEx/sampling_with_repetition/Liver/RNAseq_samples_Liver_size_100_rep_1.txt -f /path/to/GTEx/reads/rnaseq_filtered_files_by_tissue/gtex_rnaseq_Liver.gct -o /path/to/Databases/GTEx/networks/Liver/pearson_100_1.net -m pearson -n 100 -d 0.05 -p 6 -t signed -e pearson -a 0 -c bonferroni
```

This execution is done for all datasets at all sample sizes and repetitions. To facilitate the calculations, we ran them in a computational cluster using `slurm`. We used the following script to automatize the submission of the jobs:

```
scripts/gene_coexpression_networks/create_gene_coexpression_network_from_samples_cluster.py -i <input_dir> -o <output_dir> -r <rnaseq_file> -m <metric>
```


#### Consensus co-expression networks

We calculated the consensus networks between the 5 replicate genes co-expression networks from each sample size and condition. To do it, we ran the following command for each dataset:

```
Rscript scripts/gene_coexpression_networks/create_consensus_gene_coexpression_network.R -l <list_coexpression_networks_file> -n <consensus_network_file> -t <threshold> -m <method>
```

Where:

* `list_coexpression_networks_file`: file with the list of gene co-expression networks to be used in the analysis.
* `consensus_network_file`: file to save the consensus gene co-expression network.
* `threshold`: p-value threshold to consider a link as significant.
* `method`: method used to calculate the gene co-expression network.

Here we show an example of execution with the GTEx dataset for sample size 100:

```
Rscript scripts/gene_coexpression_networks/create_consensus_gene_coexpression_network.R -l /path/to/Databases/GTEx/networks/Liver/network_replicates_pearson_size_100.txt -n /path/to/Databases/GTEx/networks/Liver/consensus/consensus_pearson_100.net -t 0.05 -m pearson
```

This execution is done for all datasets at all sample sizes. To facilitate the calculations, we ran them in a computational cluster using `slurm`. We used the following script to automatize the submission of the jobs:

```
scripts/gene_coexpression_networks/create_consensus_gene_coexpression_network_cluster.py -n <networks_dir> -m <method> -t <threshold>
```


### 5. Analyze the gene co-expression networks

=> folder: `gene_coexpression_networks`

Execution:

```
Rscript scripts/gene_coexpression_networks/analyze_coexpression_network_by_significant_edges.R -c <coexpression_network_file> -o <output_results_dir> -s <output_subgraphs_dir> -f <file_name> -t <threshold> -p <ppi_file> -d <disease_genes_file> -e <essential_genes_file> -g <genes_dataset_file>
```

Where:
* `coexpression_network_file`: file with the gene co-expression network.
* `output_results_dir`: directory to save the results of the analysis.
* `output_subgraphs_dir`: directory to save the subgraphs of the gene co-expression network.
* `file_name`: name of the file to save the results of the analysis.
* `threshold`: p-value threshold to consider a link as significant.
* `ppi_file`: file with the protein-protein interaction network.
* `disease_genes_file`: file with the disease genes.
* `essential_genes_file`: file with the essential genes.
* `genes_dataset_file`: file with the genes of the gene expression dataset.

Here we show an example of execution with the GTEx dataset for sample size 100:

```
Rscript scripts/gene_coexpression_networks/analyze_coexpression_network_by_significant_edges.R -c /path/to/Databases/GTEx/networks/Liver/pearson_100_1.net -o /path/to/Databases/GTEx/analysis/Liver/ -s /path/to/Databases/GTEx/networks/Liver/subgraphs/ -f pearson_100 -t 0.05 -p /path/to/data/ppi/ppi_network.txt -d /path/to/data/disease_genes/disease_genes.txt -e /path/to/data/essential_genes/essential_genes.txt -g /path/to/Databases/GTEx/genes/gtex_genes.txt
```

This execution is done for all datasets at all sample sizes and repetitions. To facilitate the calculations, we ran them in a computational cluster using `slurm`. We used the following script to automatize the submission of the jobs:

```
scripts/gene_coexpression_networks/analyze_gene_coexpression_networks_cluster.py -i <input_dir> -p <ppi_file> -d <disease_genes_file> -e <essential_genes_file> -g <genes_dataset_file> -o <output_analysis_dir> -n <output_networks_dir>
```


### 6. Parse the results of the analyses and create summary tables

=> folder: `parse_results`

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

### 9. Analyze the gene co-expression of protein-protein interactions

==> folder: analyze_coexpression_ppi

We analyze the influence of sample size on protein-protein interaction co-expression.

The protein-protein interactions network was extracted from Gysi et al., PNAS (2023) (https://doi.org/10.1073/pnas.2301342120).

We have to run the following notebook:

```
compare_coexpression_to_ppi.Rmd
```

###  10. Create the figures

==> folder: general_analysis

We have to run the following notebook:

```
create_figures.Rmd
```
