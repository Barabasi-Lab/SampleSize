library(markdown)

# Define UI for application that draws a histogram
ui <- fluidPage(

  # Enable the use of shinyjs
  shinyjs::useShinyjs(),

  # Application title
  titlePanel("Sample size value in Network Medicine"),
  
  
  # Define navigation bar
  navbarPage("Categories",
             tabPanel("Topology",
                      
                      # Sidebar with a slider input for number of bins 
                      sidebarLayout(
                        sidebarPanel(
                          checkboxGroupInput("dataset_topology", label = "Method:", 
                                             c('GTEx' = 'gtex',
                                               'Scipher' = 'scipher',
                                               'TCGA' = 'tcga'
                                             ),
                                             selected='scipher'),
                          verbatimTextOutput(outputId = "topology_gtex_tissues"),
                          shinyWidgets::pickerInput(
                            inputId = "topology_gtex_tissues",
                            label = "Tissue (if GTEx dataset)",
                            choices = c("spleen", "whole blood"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'whole blood',
                            multiple = TRUE),
                          verbatimTextOutput(outputId = "topology_gtex_sex"),
                          shinyWidgets::pickerInput(
                            inputId = "topology_gtex_sex",
                            label = "Sex (if GTEx dataset)",
                            choices = c("female", "male", "both"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'both',
                            multiple = TRUE),
                          selectInput("type_scipher_dataset_topology", label = "Type of Scipher dataset (if Scipher dataset):", 
                                      choices = list("complete dataset" = "complete_dataset", "samples at baseline" = "samples_baseline", "sample per patient at baseline" = "sample_per_patient_baseline", "responders at baseline" = "responder_baseline", "non-responders at baseline" = "nonresponder_baseline", "old dataset (responders)" = "scipherold_responder", "old dataset (non-responders)" = "scipherold_nonresponder"), 
                                      selected = "complete_dataset"
                          ),
                        checkboxGroupInput("method_topology", label = "Method:", 
                                           c('wTO' = 'wto',
                                             'WGCNA' = 'wgcna',
                                             'ARACNE' = 'aracne',
                                             'Pearson' = 'pearson',
                                             'Spearman' = 'spearman'
                                           ), selected = 'spearman'
                        ),
                        radioButtons("parameter_topology", "Parameter:",
                                       c('Num. nodes' = 'num_nodes',
                                         'Num. edges' = 'num_edges',
                                         'Av. degree' = 'av_degree',
                                         'Av. path length' = 'av_path_length',
                                         'Av. clust. coef.' = 'av_clustering_coef',
                                         'Num. components' = 'num_components',
                                         'Nodes of the LCC' = 'num_lcc_nodes',
                                         'Edges of the LCC' = 'num_lcc_edges',
                                         'LCC z-score' = 'lcc_z',
                                         'LCC p-value' = 'lcc_pvalue',
                                         'LCC log(p-value)' = 'log_lcc_pvalue',
                                         'Maximum k-core' = 'max_k',
                                         'Num. nodes in main core' = 'num_main_core_nodes',
                                         'Num. edges in main core' = 'num_main_core_edges',
                                         'Repetitions' = 'rep')),
                        checkboxGroupInput("pvalue_threshold_topology", label = "P-value threshold:", 
                                           choices = c('0.05' = "0.05-NA",
                                                       '0.001' = "0.001-NA",
                                                       '0.05, disparity 0.05' = "0.05-0.05",
                                                       '0.001, disparity 0.05' = "0.001-0.05"),
                                           selected = "0.05-NA",
                                           ),
                          shinyWidgets::materialSwitch(inputId = "log_x_topology", label = "Turn x axis to log scale:", status="primary"),
                          shinyWidgets::materialSwitch(inputId = "log_y_topology", label = "Turn y axis to log scale:", status="primary"),
                          shinyWidgets::materialSwitch(inputId = "sd_topology", label = "Plot standard deviation:", status="primary")
                        ),
                        
                        # Show boxplot of topological parameters
                        mainPanel(
                          plotOutput("topologyBoxPlot")
                        )
                      )
             ),

             
             tabPanel("PPI",
                      
                      # Sidebar with a slider input for number of bins 
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons("dataset_ppi", "Dataset:",
                                       c('GTEx' = 'gtex',
                                         'Scipher' = 'scipher',
                                         'TCGA' = 'tcga'
                                       ),
                                       selected='scipher'),
                          verbatimTextOutput(outputId = "ppi_gtex_tissues"),
                          shinyWidgets::pickerInput(
                            inputId = "ppi_gtex_tissues",
                            label = "Tissue (if GTEx dataset)",
                            choices = c("spleen", "whole blood"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'whole blood',
                            multiple = TRUE),
                          verbatimTextOutput(outputId = "ppi_gtex_sex"),
                          shinyWidgets::pickerInput(
                            inputId = "ppi_gtex_sex",
                            label = "Sex (if GTEx dataset)",
                            choices = c("female", "male", "both"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'both',
                            multiple = TRUE),
                          selectInput("type_scipher_dataset_ppi", label = "Type of Scipher dataset (if Scipher dataset):", 
                                      choices = list("complete dataset" = "complete_dataset", "samples at baseline" = "samples_baseline", "sample per patient at baseline" = "sample_per_patient_baseline", "responders at baseline" = "responder_baseline", "non-responders at baseline" = "nonresponder_baseline", "old dataset (responders)" = "scipherold_responder", "old dataset (non-responders)" = "scipherold_nonresponder"), 
                                      selected = "complete_dataset"
                          ),
                          checkboxGroupInput("method_ppi", label = "Method:", 
                                             c('wTO' = 'wto',
                                               'WGCNA' = 'wgcna',
                                               'ARACNE' = 'aracne',
                                               'Pearson' = 'pearson',
                                               'Spearman' = 'spearman'
                                             ), selected = 'spearman'
                          ),
                          radioButtons("parameter_ppi", "Parameter:",
                                       c('Num. PPI nodes' = 'num_ppi_nodes',
                                         'Num. PPI edges' = 'num_ppi_edges',
                                         'Fraction PPI nodes' = 'fraction_ppi_nodes',
                                         'Fraction PPI edges' = 'fraction_ppi_edges',
                                         'Num. PPI main core nodes' = 'num_ppi_main_core_nodes',
                                         'Num. PPI main core edges' = 'num_ppi_main_core_edges',
                                         'Fraction PPI main core nodes' = 'fraction_ppi_main_core_nodes',
                                         'Fraction PPI main core edges' = 'fraction_ppi_main_core_edges')),
                          checkboxGroupInput("pvalue_threshold_ppi", label = "P-value threshold:", 
                                             choices = c('0.05' = "0.05-NA",
                                                         '0.001' = "0.001-NA",
                                                         '0.05, disparity 0.05' = "0.05-0.05",
                                                         '0.001, disparity 0.05' = "0.001-0.05"),
                                             selected = "0.05-NA",
                          ),
                          shinyWidgets::materialSwitch(inputId = "log_x_ppi", label = "Turn x axis to log scale:", status="primary"),
                          shinyWidgets::materialSwitch(inputId = "log_y_ppi", label = "Turn y axis to log scale:", status="primary")
                        ),
                        
                        # Show boxplot of topological parameters
                        mainPanel(
                          plotOutput("ppiBoxPlot")
                        )
                      )
             ),
             
             
             tabPanel("Essentiality",
                      
                      # Sidebar with a slider input for number of bins 
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons("dataset_essentiality", "Dataset:",
                                       c('GTEx' = 'gtex',
                                         'Scipher' = 'scipher',
                                         'TCGA' = 'tcga'
                                       ),
                                       selected='scipher'),
                          verbatimTextOutput(outputId = "essentiality_gtex_tissues"),
                          shinyWidgets::pickerInput(
                            inputId = "essentiality_gtex_tissues",
                            label = "Tissue (if GTEx dataset)",
                            choices = c("spleen", "whole blood"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'whole blood',
                            multiple = TRUE),
                          verbatimTextOutput(outputId = "essentiality_gtex_sex"),
                          shinyWidgets::pickerInput(
                            inputId = "essentiality_gtex_sex",
                            label = "Sex (if GTEx dataset)",
                            choices = c("female", "male", "both"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'both',
                            multiple = TRUE),
                          selectInput("type_scipher_dataset_essentiality", label = "Type of Scipher dataset (if Scipher dataset):", 
                                      choices = list("complete dataset" = "complete_dataset", "samples at baseline" = "samples_baseline", "sample per patient at baseline" = "sample_per_patient_baseline", "responders at baseline" = "responder_baseline", "non-responders at baseline" = "nonresponder_baseline", "old dataset (responders)" = "scipherold_responder", "old dataset (non-responders)" = "scipherold_nonresponder"), 
                                      selected = "all_samples"
                          ),
                          radioButtons("method_essentiality", "Method:",
                                       c('wTO' = 'wto',
                                         'WGCNA' = 'wgcna',
                                         'ARACNE' = 'aracne',
                                         'Pearson' = 'pearson',
                                         'Spearman' = 'spearman'
                                       ), selected = 'spearman'
                                       ),
                          radioButtons("parameter_essentiality", "Parameter:",
                                       c('Num. essential nodes' = 'num_essential_genes',
                                         'Fraction of essential genes' = 'fraction_essential_genes',
                                         'Num. components' = 'num_components',
                                         'Num. nodes of the LCC' = 'num_essential_lcc_nodes',
                                         'Essential genes rLCC' = 'fraction_essential_lcc_nodes',
                                         'Num. edges of the LCC' = 'num_lcc_edges',
                                         'LCC z-score' = 'lcc_z',
                                         'LCC p-value' = 'lcc_pvalue',
                                         'LCC log(p-value)' = 'log_lcc_pvalue')),
                          radioButtons("pvalue_threshold_essentiality", "P-value threshold",
                                       c('0.05' = "0.05-NA",
                                         '0.001' = "0.001-NA",
                                         '0.05, disparity 0.05' = "0.05-0.05",
                                         '0.001, disparity 0.05' = "0.001-0.05",
                                         'Group' = 'group_essentiality_pvalue_threshold')),
                          shinyWidgets::materialSwitch(inputId = "log_essentiality", label = "Turn log scale:", status="primary")
                        ),
                        
                        # Show boxplot of topological parameters
                        mainPanel(
                          plotOutput("essentialityBoxPlot")
                        )
                      )
             ),
             
             
             navbarMenu("Disease genes",
                        
                        tabPanel("By thresholds",
                                 
                                 # Sidebar with a slider input for number of bins 
                                 sidebarLayout(
                                   sidebarPanel(
                                     radioButtons("dataset_disease_genes", "Dataset:",
                                                  c('GTEx' = 'gtex',
                                                    'Scipher' = 'scipher',
                                                    'TCGA' = 'tcga'
                                                  ), selected = "scipher"
                                                  ),
                                     shinyWidgets::pickerInput(
                                       inputId = "disease_genes_gtex_tissues",
                                       label = "Tissue (if GTEx dataset)",
                                       choices = c("spleen", "whole blood"),
                                       options = list(
                                         `actions-box` = TRUE,
                                         size = 10,
                                         `selected-text-format` = "count > 1"
                                       ),
                                       selected = 'whole blood',
                                       multiple = TRUE
                                     ),
                                     verbatimTextOutput(outputId = "disease_genes_gtex_tissues"),
                                     shinyWidgets::pickerInput(
                                       inputId = "disease_genes_gtex_sex",
                                       label = "Sex (if GTEx dataset)",
                                       choices = c("female", "male", "both"),
                                       options = list(
                                         `actions-box` = TRUE,
                                         size = 10,
                                         `selected-text-format` = "count > 1"
                                       ),
                                       selected = 'both',
                                       multiple = TRUE
                                     ),
                                     verbatimTextOutput(outputId = "disease_genes_gtex_sex"),
                                     selectInput("type_scipher_dataset_disease_genes", label = "Type of Scipher dataset (if Scipher dataset):", 
                                                 choices = list("complete dataset" = "complete_dataset", "samples at baseline" = "samples_baseline", "sample per patient at baseline" = "sample_per_patient_baseline", "responders at baseline" = "responder_baseline", "non-responders at baseline" = "nonresponder_baseline", "old dataset (responders)" = "scipherold_responder", "old dataset (non-responders)" = "scipherold_nonresponder"), 
                                                 selected = "all_samples"
                                     ),
                                     radioButtons("method_disease_genes", "Method:",
                                                  c('wTO' = 'wto',
                                                    'WGCNA' = 'wgcna',
                                                    'ARACNE' = 'aracne',
                                                    'Pearson' = 'pearson',
                                                    'Spearman' = 'spearman'
                                                  ), selected = 'spearman'
                                                  ),
                                     radioButtons("parameter_disease_genes", "Parameter:",
                                                  c('Num. disease genes in network' = 'num_disease_genes',
                                                    'Fraction disease genes in network' = 'fraction_disease_genes',
                                                    'Num. components of disease genes' = 'num_disease_components',
                                                    'Num. disease genes forming a LCC' = 'num_disease_lcc_nodes',
                                                    'Disease rLCC' = 'fraction_disease_lcc_nodes',
                                                    'Num. edges of the disease gene LCC' = 'num_disease_lcc_edges',
                                                    'Disease gene LCC z-score' = 'disease_lcc_z',
                                                    'Disease gene LCC p-value' = 'disease_lcc_pvalue',
                                                    'Disease gene LCC log(p-value)' = 'log_disease_lcc_pvalue')),
                                 radioButtons("pvalue_threshold_disease_genes", "P-value threshold",
                                              c('0.05' = "0.05-NA",
                                                '0.001' = "0.001-NA",
                                                '0.05, disparity 0.05' = "0.05-0.05",
                                                '0.001, disparity 0.05' = "0.001-0.05")),
                                     radioButtons("group_by", "Group by:",
                                                  c('Disease' = 'group_disease',
                                                    'Disease class' = 'group_disease_class',
                                                    'P-value' = 'group_disease_pvalue',
                                                    'Not group' = 'no_group')),
                                     shinyWidgets::pickerInput(
                                       inputId = "immune",
                                       label = "Immune system diseases",
                                       choices = c("immune system diseases", "arthritis rheumatoid", "asthma", "graves disease", "multiple sclerosis"),
                                       options = list(
                                         `actions-box` = TRUE,
                                         size = 10,
                                         `selected-text-format` = "count > 1"
                                       ),
                                       selected = 'arthritis rheumatoid',
                                       multiple = TRUE
                                     ),
                                     verbatimTextOutput(outputId = "diseases_immune"),
                                     shinyWidgets::pickerInput(
                                       inputId = "cardiovascular",
                                       label = "Cardiovascular diseases",
                                       choices = c("cardiovascular diseases", "arrhythmias, cardiac", "cardiomyopathies", "heart arrest", "myocardial ischemia"),
                                       options = list(
                                         `actions-box` = TRUE,
                                         size = 10,
                                         `selected-text-format` = "count > 1"
                                       ),
                                       #selected = 'heart arrest',
                                       multiple = TRUE
                                     ),
                                     verbatimTextOutput(outputId = "diseases_cardiovascular"),
                                     shinyWidgets::pickerInput(
                                       inputId = "nervous",
                                       label = "Nervous system diseases",
                                       choices = c("nervous system diseases", "alzheimer disease", "epilepsy", "parkinson disease"),
                                       options = list(
                                         `actions-box` = TRUE,
                                         size = 10,
                                         `selected-text-format` = "count > 1"
                                       ),
                                       #selected = 'alzheimer disease',
                                       multiple = TRUE
                                     ),
                                     verbatimTextOutput(outputId = "diseases_nervous"),
                                     shinyWidgets::pickerInput(
                                       inputId = "neoplasms",
                                       label = "Neoplasms",
                                       choices = c("neoplasms", "breast neoplasms", "colorectal neoplasm", "kidney neoplasms", "leukemia, myeloid", "lung neoplasms"),
                                       options = list(
                                         `actions-box` = TRUE,
                                         size = 10,
                                         `selected-text-format` = "count > 1"
                                       ),
                                       #selected = 'breast neoplasms',
                                       multiple = TRUE
                                     ),
                                     verbatimTextOutput(outputId = "diseases_neoplasms"),
                                     shinyWidgets::pickerInput(
                                       inputId = "digestive",
                                       label = "Digestive system diseases",
                                       choices = c("digestive system diseases", "crohn disease", "cholestasis", "liver cirrhosis", "gastroenteritis"),
                                       options = list(
                                         `actions-box` = TRUE,
                                         size = 10,
                                         `selected-text-format` = "count > 1"
                                       ),
                                       #selected = 'liver cirrhosis',
                                       multiple = TRUE
                                     ),
                                     verbatimTextOutput(outputId = "diseases_digestive"),
                                     shinyWidgets::pickerInput(
                                       inputId = "nutritional",
                                       label = "Nutritional and metabolic diseases",
                                       choices = c("nutritional and metabolic diseases", "amyotrophic lateral sclerosis", "celiac disease", "diabetes mellitus type 1", "diabetes mellitus type 2", "obesity"),
                                       options = list(
                                         `actions-box` = TRUE,
                                         size = 10,
                                         `selected-text-format` = "count > 1"
                                       ),
                                       #selected = 'diabetes mellitus type 2',
                                       multiple = TRUE
                                     ),
                                     verbatimTextOutput(outputId = "diseases_nutritional"),
                                     shinyWidgets::pickerInput(
                                       inputId = "endocrine",
                                       label = "Endocrine system diseases",
                                       choices = c("endocrine system diseases", "diabetes mellitus type 2", "hyperthyroidism", "goiter"),
                                       options = list(
                                         `actions-box` = TRUE,
                                         size = 10,
                                         `selected-text-format` = "count > 1"
                                       ),
                                       #selected = 'hyperthyroidism',
                                       multiple = TRUE
                                     ),
                                     verbatimTextOutput(outputId = "diseases_endocrine"),
                                     shinyWidgets::materialSwitch(inputId = "log_x_diseases", label = "Turn x axis to log scale:", status="primary"),
                                     shinyWidgets::materialSwitch(inputId = "log_y_diseases", label = "Turn y axis to log scale:", status="primary")
                                   ),
                                   
                                   # Show boxplot of topological parameters
                                   mainPanel(
                                     plotOutput("diseaseBoxPlot")
                                   )
                                 )
                        ),
                        
                        tabPanel("By score",
                                 
                                 # Sidebar with a slider input for number of bins 
                                 sidebarLayout(
                                   sidebarPanel(
                                     radioButtons("dataset_disease_genes_by_score", "Dataset:",
                                                  c('GTEx' = 'gtex',
                                                    'Scipher' = 'scipher',
                                                    'TCGA' = 'tcga'
                                                  )),
                                     radioButtons("method_disease_genes_by_score", "Method:",
                                                  c('wTO' = 'wto',
                                                    'WGCNA' = 'wgcna',
                                                    'ARACNE' = 'aracne',
                                                    'Pearson' = 'pearson',
                                                    'Spearman' = 'spearman'
                                                  )),
                                     radioButtons("metric_disease_genes_by_score", "Evaluation metric:",
                                                  c('Difference' = 'difference',
                                                    'Number of edges' = 'num.edges',
                                                    'Co-expression (network from subset of samples)' = 'score.subset.mean',
                                                    'Co-expression (network from all samples)' = 'score.all.samples.mean',
                                                    'Distance' = 'distances',
                                                    'Local clustering coeff.' = 'local_clustering_coefficient.mean',
                                                    'Degree' = 'degree_centrality.mean',
                                                    'Betweenness centrality' = 'betweenness_centrality.mean'
                                                  )),
                                     radioButtons("disease_genes_by_score_group_by", "Group by:",
                                                  c('Disease' = 'group_disease',
                                                    'Disease class' = 'group_disease_class',
                                                    'P-value' = 'group_disease_pval',
                                                    '% of top scoring edges' = 'group_disease_top',
                                                    'Not group' = 'no_group')),
                                     shinyWidgets::pickerInput(
                                       inputId = "disease_genes_by_score_immune",
                                       label = "Immune system diseases",
                                       choices = c("immune system diseases", "arthritis, rheumatoid", "asthma", "graves disease", "multiple sclerosis"),
                                       options = list(
                                         `actions-box` = TRUE,
                                         size = 10,
                                         `selected-text-format` = "count > 1"
                                       ),
                                       selected = 'arthritis, rheumatoid',
                                       multiple = TRUE
                                     ),
                                     verbatimTextOutput(outputId = "disease_genes_by_score_immune"),
                                     shinyWidgets::pickerInput(
                                       inputId = "disease_genes_by_score_cardiovascular",
                                       label = "Cardiovascular diseases",
                                       choices = c("cardiovascular diseases", "arrhythmias, cardiac", "cardiomyopathies", "heart arrest", "myocardial ischemia"),
                                       options = list(
                                         `actions-box` = TRUE,
                                         size = 10,
                                         `selected-text-format` = "count > 1"
                                       ),
                                       selected = 'arrhythmias, cardiac',
                                       multiple = TRUE
                                     ),
                                     verbatimTextOutput(outputId = "disease_genes_by_score_cardiovascular"),
                                     shinyWidgets::pickerInput(
                                       inputId = "disease_genes_by_score_nervous",
                                       label = "Nervous system diseases",
                                       choices = c("nervous system diseases", "alzheimer disease", "epilepsy", "parkinson disease"),
                                       options = list(
                                         `actions-box` = TRUE,
                                         size = 10,
                                         `selected-text-format` = "count > 1"
                                       ),
                                       selected = 'alzheimer disease',
                                       multiple = TRUE
                                     ),
                                     verbatimTextOutput(outputId = "disease_genes_by_score_nervous"),
                                     shinyWidgets::pickerInput(
                                       inputId = "disease_genes_by_score_neoplasms",
                                       label = "Neoplasms",
                                       choices = c("neoplasms", "breast neoplasms", "colorectal neoplasm", "kidney neoplasms", "leukemia, myeloid", "lung neoplasms"),
                                       options = list(
                                         `actions-box` = TRUE,
                                         size = 10,
                                         `selected-text-format` = "count > 1"
                                       ),
                                       selected = 'breast neoplasms',
                                       multiple = TRUE
                                     ),
                                     verbatimTextOutput(outputId = "disease_genes_by_score_neoplasms"),
                                     shinyWidgets::pickerInput(
                                       inputId = "disease_genes_by_score_digestive",
                                       label = "Digestive system diseases",
                                       choices = c("digestive system diseases", "crohn disease", "cholestasis", "liver cirrhosis", "gastroenteritis"),
                                       options = list(
                                         `actions-box` = TRUE,
                                         size = 10,
                                         `selected-text-format` = "count > 1"
                                       ),
                                       #selected = 'liver cirrhosis',
                                       multiple = TRUE
                                     ),
                                     verbatimTextOutput(outputId = "disease_genes_by_score_digestive"),
                                     shinyWidgets::pickerInput(
                                       inputId = "disease_genes_by_score_nutritional",
                                       label = "Nutritional and metabolic diseases",
                                       choices = c("nutritional and metabolic diseases", "amyotrophic lateral sclerosis", "celiac disease", "diabetes mellitus type 1", "diabetes mellitus type 2", "obesity"),
                                       options = list(
                                         `actions-box` = TRUE,
                                         size = 10,
                                         `selected-text-format` = "count > 1"
                                       ),
                                       #selected = 'diabetes mellitus type 2',
                                       multiple = TRUE
                                     ),
                                     verbatimTextOutput(outputId = "disease_genes_by_score_nutritional"),
                                     shinyWidgets::pickerInput(
                                       inputId = "disease_genes_by_score_endocrine",
                                       label = "Endocrine system diseases",
                                       choices = c("endocrine system diseases", "diabetes mellitus type 2", "hyperthyroidism", "goiter"),
                                       options = list(
                                         `actions-box` = TRUE,
                                         size = 10,
                                         `selected-text-format` = "count > 1"
                                       ),
                                       #selected = 'hyperthyroidism',
                                       multiple = TRUE
                                     ),
                                     verbatimTextOutput(outputId = "disease_genes_by_score_endocrine"),
                                     
                                     shinyWidgets::pickerInput(
                                       inputId = "disease_genes_by_score_gtex_tissues",
                                       label = "Tissue (if GTEx dataset)",
                                       choices = c("spleen", "whole blood"),
                                       options = list(
                                         `actions-box` = TRUE,
                                         size = 10,
                                         `selected-text-format` = "count > 1"
                                       ),
                                       selected = 'whole blood',
                                       multiple = TRUE
                                     ),
                                     verbatimTextOutput(outputId = "disease_genes_by_score_gtex_tissues"),
                                     shinyWidgets::pickerInput(
                                       inputId = "disease_genes_by_score_gtex_sex",
                                       label = "Sex (if GTEx dataset)",
                                       choices = c("female", "male", "both"),
                                       options = list(
                                         `actions-box` = TRUE,
                                         size = 10,
                                         `selected-text-format` = "count > 1"
                                       ),
                                       selected = 'both',
                                       multiple = TRUE
                                     ),
                                     verbatimTextOutput(outputId = "disease_genes_by_score_gtex_sex"),
                                     selectInput("type_scipher_dataset_disease_genes_by_score", label = "Type of Scipher dataset (if Scipher dataset):", 
                                                 choices = list("complete dataset" = "complete_dataset", "samples at baseline" = "samples_baseline", "sample per patient at baseline" = "sample_per_patient_baseline", "responders at baseline" = "responder_baseline", "non-responders at baseline" = "nonresponder_baseline", "old dataset (responders)" = "scipherold_responder", "old dataset (non-responders)" = "scipherold_nonresponder"), 
                                                 selected = "all_samples"
                                     )
                                     
                                   ),
                                   
                                   # Show boxplot of topological parameters
                                   mainPanel(
                                     plotOutput("diseaseByScoreBoxPlot")
                                   )
                                 )
                        )
             ),
             
             
             
             tabPanel("Similarity",
                      
                      # Sidebar with a slider input for number of bins 
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons("dataset_similarity", "Dataset:",
                                       c('GTEx' = 'gtex',
                                         'Scipher' = 'scipher',
                                         'TCGA' = 'tcga'
                                       ),
                                       selected='scipher'),
                          verbatimTextOutput(outputId = "similarity_gtex_tissues"),
                          shinyWidgets::pickerInput(
                            inputId = "tissue",
                            label = "Tissue (if GTEx dataset)",
                            choices = c("spleen", "whole blood"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'whole blood',
                            multiple = TRUE),
                          verbatimTextOutput(outputId = "similarity_gtex_sex"),
                          shinyWidgets::pickerInput(
                            inputId = "sex",
                            label = "Sex (if GTEx dataset)",
                            choices = c("female", "male", "both"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'both',
                            multiple = TRUE),
                          selectInput("type_scipher_dataset_similarity", label = "Type of Scipher dataset (if Scipher dataset):", 
                                      choices = list("complete dataset" = "complete_dataset", "samples at baseline" = "samples_baseline", "sample per patient at baseline" = "sample_per_patient_baseline", "responders at baseline" = "responder_baseline", "non-responders at baseline" = "nonresponder_baseline", "old dataset (responders)" = "scipherold_responder", "old dataset (non-responders)" = "scipherold_nonresponder"), 
                                      selected = "all_samples"
                          ),
                          radioButtons("method_similarity", "Method:",
                                       c('wTO' = 'wto',
                                         'WGCNA' = 'wgcna',
                                         'ARACNE' = 'aracne',
                                         'Pearson' = 'pearson',
                                         'Spearman' = 'spearman'
                                       ), selected = 'spearman'
                          ),
                          radioButtons("parameter_similarity", "Parameter:",
                                       c('Overlap index' = 'overlapindex',
                                         'Jaccard index' = 'jaccardIndex'
                                         ), selected = 'overlapindex'),
                          radioButtons("type_analysis_similarity", "Type of analysis:",
                                       c('Network' = 'subgraphs',
                                         'Main core' = 'main_core',
                                         'Essential genes' = 'essential_genes',
                                         'RA genes' = 'disease_genes'
                                       ), selected = 'subgraphs'),
                          radioButtons("pvalue_threshold_similarity", "P-value threshold",
                                       c('0.05' = 0.05,
                                         '0.001' = 0.001))
                        ),
                        
                        # Show boxplot of topological parameters
                        fluidRow(
                          column(6, plotOutput("similarityBoxPlot")),
                          column(6, plotOutput("comparisonBoxPlot"))
                        )
                      )
             ),
             
             
             
             tabPanel("Co-expressed PPIs",
                      
                      # Sidebar with a slider input for number of bins 
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons("dataset_coexpressed_ppis", "Dataset:",
                                       c('GTEx' = 'gtex',
                                         'Scipher' = 'scipher',
                                         'TCGA' = 'tcga'
                                       ), selected = 'scipher'
                                       ),
                          shinyWidgets::pickerInput(
                            inputId = "coexpressed_ppis_gtex_tissues",
                            label = "Tissue (if GTEx dataset)",
                            choices = c("spleen", "whole blood"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'whole blood',
                            multiple = TRUE
                          ),
                          verbatimTextOutput(outputId = "coexpressed_ppis_gtex_tissues"),
                          shinyWidgets::pickerInput(
                            inputId = "coexpressed_ppis_gtex_sex",
                            label = "Sex (if GTEx dataset)",
                            choices = c("female", "male", "both"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'both',
                            multiple = TRUE
                          ),
                          verbatimTextOutput(outputId = "coexpressed_ppis_gtex_sex"),
                          selectInput("type_scipher_dataset_coexpressed_ppis", label = "Type of Scipher dataset (if Scipher dataset):", 
                                      choices = list("complete dataset" = "complete_dataset", "samples at baseline" = "samples_baseline", "sample per patient at baseline" = "sample_per_patient_baseline", "responders at baseline" = "responder_baseline", "non-responders at baseline" = "nonresponder_baseline", "old dataset (responders)" = "scipherold_responder", "old dataset (non-responders)" = "scipherold_nonresponder"), 
                                      selected = "all_samples"
                          ),
                          radioButtons("method_coexpressed_ppis", "Method:",
                                       c('wTO' = 'wto',
                                         'WGCNA' = 'wgcna',
                                         'ARACNE' = 'aracne',
                                         'Pearson' = 'pearson',
                                         'Spearman' = 'spearman'
                                       ), selected = 'wgcna'
                                       ),
                          radioButtons("metric_coexpressed_ppis", "Evaluation metric:",
                                       c('Accuracy' = 'accuracy',
                                         'F1' = 'F1',
                                         'MCC' = 'MCC',
                                         'Correlation' = 'corr',
                                         'TPR' = 'TPR',
                                         'FPR' = 'FPR',
                                         'TP' = 'TP',
                                         'FP' = 'FP',
                                         'TN' = 'TN',
                                         'FN' = 'FN')),
                          radioButtons("parameter_coexpressed_ppis", "Parameter to predict high co-expression:",
                                       c('Short distance' = 'balanced',
                                         'Short distance & high clustering' = 'clustering',
                                         'Short distance & high degree' = 'degree',
                                         'Short distance & high betweenness' = 'betweenness'))
                        ),
                          
                        # Show boxplot of topological parameters
                        mainPanel(
                          plotOutput("coexpressedPPIsBoxPlot")
                        )
                      )
             ),


             tabPanel("Stability",
                      
                      # Sidebar with a slider input for number of bins 
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons("dataset_stability", "Dataset:",
                                       c('GTEx' = 'gtex',
                                         'Scipher' = 'scipher',
                                         'TCGA' = 'tcga'
                                       )),
                          shinyWidgets::pickerInput(
                            inputId = "stability_gtex_tissues",
                            label = "Tissue (if GTEx dataset)",
                            choices = c("spleen", "whole blood"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'whole blood',
                            multiple = TRUE
                          ),
                          verbatimTextOutput(outputId = "stability_gtex_tissues"),
                          shinyWidgets::pickerInput(
                            inputId = "stability_gtex_sex",
                            label = "Sex (if GTEx dataset)",
                            choices = c("female", "male", "both"),
                            options = list(
                              `actions-box` = TRUE,
                              size = 10,
                              `selected-text-format` = "count > 1"
                            ),
                            selected = 'both',
                            multiple = TRUE
                          ),
                          verbatimTextOutput(outputId = "stability_gtex_sex"),
                          selectInput("type_scipher_dataset_stability", label = "Type of Scipher dataset (if Scipher dataset):", 
                                      choices = list("complete dataset" = "complete_dataset", "samples at baseline" = "samples_baseline", "sample per patient at baseline" = "sample_per_patient_baseline", "responders at baseline" = "responder_baseline", "non-responders at baseline" = "nonresponder_baseline", "old dataset (responders)" = "scipherold_responder", "old dataset (non-responders)" = "scipherold_nonresponder"), 
                                      selected = "all_samples"
                          ),
                          radioButtons("method_stability", "Method:",
                                       c('wTO' = 'wto',
                                         'WGCNA' = 'wgcna',
                                         'ARACNE' = 'aracne',
                                         'Pearson' = 'pearson',
                                         'Spearman' = 'spearman'
                                       )),
                          radioButtons("type_interaction_stability", "Types of protein pairs (from the PPI network) to analyze:",
                                       c('All protein pairs' = 'all',
                                         'Only direct interactions' = 'direct',
                                         'Only pairs not interacting directly' = 'indirect')),
                          radioButtons("metric_stability", "Evaluation metric:",
                                       c('Difference' = 'difference.mean',
                                         'Number of edges' = 'number.of.edges',
                                         'Co-expression (network from subset of samples)' = 'score.subset.mean',
                                         'Co-expression (network from all samples)' = 'score.all.samples.mean',
                                         'Distance' = 'distance.mean',
                                         'Local clustering coeff.' = 'local_clustering_coefficient.mean',
                                         'Degree' = 'degree_centrality.mean',
                                         'Betweenness centrality' = 'betweenness_centrality.mean'
                                         )),
                          radioButtons("difference_range_stability", "Co-expression score difference:",
                                       c('> 0' = '0-Inf',
                                         '0 - 0.001' = '0-0.001',
                                         '0.001 - 0.05' = '0.001-0.05',
                                         '0.05 - 0.1' = '0.05-0.1',
                                         '> 0.1' = '0.1-Inf'
                                         )),
                          radioButtons("stability_group_by", "Group by:",
                                       c('Type of interaction' = 'group_type_interaction',
                                         'Range of co-expression score difference' = 'group_difference_range',
                                         'Not group' = 'no_group'))
                        ),
                        
                        # Show boxplot of topological parameters
                        mainPanel(
                          plotOutput("stabilityBoxPlot")
                        )
                      )
             ),

             navbarMenu("Co-expression scores",
                        
               tabPanel("By ranges",
                        
                        # Sidebar with a slider input for number of bins 
                        sidebarLayout(
                          sidebarPanel(
                            radioButtons("dataset_scores", "Dataset:",
                                         c('GTEx' = 'gtex',
                                           'Scipher' = 'scipher',
                                           'TCGA' = 'tcga'
                                         )),
                            shinyWidgets::pickerInput(
                              inputId = "scores_gtex_tissues",
                              label = "Tissue (if GTEx dataset)",
                              choices = c("spleen", "whole blood"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              selected = 'whole blood',
                              multiple = TRUE
                            ),
                            verbatimTextOutput(outputId = "scores_gtex_tissues"),
                            shinyWidgets::pickerInput(
                              inputId = "scores_gtex_sex",
                              label = "Sex (if GTEx dataset)",
                              choices = c("female", "male", "both"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              selected = 'both',
                              multiple = TRUE
                            ),
                            verbatimTextOutput(outputId = "scores_gtex_sex"),
                            selectInput("type_scipher_dataset_scores", label = "Type of Scipher dataset (if Scipher dataset):", 
                                        choices = list("complete dataset" = "complete_dataset", "samples at baseline" = "samples_baseline", "sample per patient at baseline" = "sample_per_patient_baseline", "responders at baseline" = "responder_baseline", "non-responders at baseline" = "nonresponder_baseline", "old dataset (responders)" = "scipherold_responder", "old dataset (non-responders)" = "scipherold_nonresponder"), 
                                        selected = "all_samples"
                            ),
                            radioButtons("method_scores", "Method:",
                                         c('wTO' = 'wto',
                                           'WGCNA' = 'wgcna',
                                           'ARACNE' = 'aracne',
                                           'Pearson' = 'pearson',
                                           'Spearman' = 'spearman'
                                         )),
                            radioButtons("type_interaction_scores", "Types of protein pairs (from the PPI network) to analyze:",
                                         c('All protein pairs' = 'all',
                                           'Only direct interactions' = 'direct',
                                           'Only pairs not interacting directly' = 'indirect',
                                           'Group by type of interaction' = 'group')),
                            radioButtons("metric_scores", "Evaluation metric:",
                                         c('Difference' = 'difference.mean',
                                           'Number of edges' = 'number.of.edges',
                                           'Co-expression (network from subset of samples)' = 'score.subset.mean',
                                           'Co-expression (network from all samples)' = 'score.all.samples.mean',
                                           'Distance' = 'distance.mean',
                                           'Local clustering coeff.' = 'local_clustering_coefficient.mean',
                                           'Degree' = 'degree_centrality.mean',
                                           'Betweenness centrality' = 'betweenness_centrality.mean'
                                         )),
                            checkboxGroupInput("score_range", label = "Co-expression difference:", 
                                               choices = list("All ranges grouped" = "group", 
                                                              "-1 to -0.9" = "-1--0.9", 
                                                              "-0.9 to -0.8" = "-0.9--0.8",
                                                              "0 to 0.1" = "0-0.1",
                                                              "0.1 to 0.2" = "0.1-0.2",
                                                              "0.2 to 0.3" = "0.2-0.3",
                                                              "0.3 to 0.4" = "0.3-0.4",
                                                              "0.4 to 0.5" = "0.4-0.5",
                                                              "0.5 to 0.6" = "0.5-0.6",
                                                              "0.6 to 0.7" = "0.6-0.7",
                                                              "0.7 to 0.8" = "0.7-0.8",
                                                              "0.8 to 0.9" = "0.8-0.9",
                                                              "0.9 to 1" = "0.9-1"
                                               ),
                                               selected = "group")
                          ),
                          
                          # Show boxplot of topological parameters
                          mainPanel(
                            plotOutput("scoresBoxPlot")
                          )
                        )
               ),
               
               tabPanel("By thresholds",
                        
                        # Sidebar with a slider input for number of bins 
                        sidebarLayout(
                          sidebarPanel(
                            radioButtons("dataset_scores_thresholds", "Dataset:",
                                         c('GTEx' = 'gtex',
                                           'Scipher' = 'scipher',
                                           'TCGA' = 'tcga'
                                         )),
                            shinyWidgets::pickerInput(
                              inputId = "scores_thresholds_gtex_tissues",
                              label = "Tissue (if GTEx dataset)",
                              choices = c("spleen", "whole blood"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              selected = 'whole blood',
                              multiple = TRUE
                            ),
                            verbatimTextOutput(outputId = "scores_thresholds_gtex_tissues"),
                            shinyWidgets::pickerInput(
                              inputId = "scores_thresholds_gtex_sex",
                              label = "Sex (if GTEx dataset)",
                              choices = c("female", "male", "both"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              selected = 'both',
                              multiple = TRUE
                            ),
                            verbatimTextOutput(outputId = "scores_thresholds_gtex_sex"),
                            selectInput("type_scipher_dataset_scores_thresholds", label = "Type of Scipher dataset (if Scipher dataset):", 
                                        choices = list("complete dataset" = "complete_dataset", "samples at baseline" = "samples_baseline", "sample per patient at baseline" = "sample_per_patient_baseline", "responders at baseline" = "responder_baseline", "non-responders at baseline" = "nonresponder_baseline", "old dataset (responders)" = "scipherold_responder", "old dataset (non-responders)" = "scipherold_nonresponder"), 
                                        selected = "all_samples"
                            ),
                            radioButtons("method_scores_thresholds", "Method:",
                                         c('wTO' = 'wto',
                                           'WGCNA' = 'wgcna',
                                           'ARACNE' = 'aracne',
                                           'Pearson' = 'pearson',
                                           'Spearman' = 'spearman'
                                         )),
                            radioButtons("type_interaction_scores_thresholds", "Types of protein pairs (from the PPI network) to analyze:",
                                         c('All protein pairs' = 'all',
                                           'Only direct interactions' = 'direct',
                                           'Only pairs not interacting directly' = 'indirect')),
                            radioButtons("metric_scores_thresholds", "Evaluation metric:",
                                         c('Number of edges (network from subset of samples)' = 'number.of.edges.subset',
                                           'Number of edges (network from all samples)' = 'number.of.edges.all.samples',
                                           'Number of edges overlapped' = 'number.of.edges.overlapped',
                                           'Number of edges lost (from the network of all samples)' = 'number.of.edges.lost',
                                           'Number of edges gained (new in the subset network)' = 'number.of.edges.gained'
                                         )),
                            checkboxGroupInput("score_thresholds", label = "Score threshold:", 
                                               choices = list("All ranges grouped" = "group", 
                                                              "0.1" = "0.1",
                                                              "0.2" = "0.2",
                                                              "0.3" = "0.3",
                                                              "0.4" = "0.4",
                                                              "0.5" = "0.5",
                                                              "0.6" = "0.6",
                                                              "0.7" = "0.7",
                                                              "0.8" = "0.8",
                                                              "0.9" = "0.9"
                                               ),
                                               selected = "group")
                          ),
                          
                          # Show boxplot of topological parameters
                          mainPanel(
                            plotOutput("scoresThresholdsBoxPlot")
                          )
                        )
               )
             )
             
             
  )
  
)
