library(markdown)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Enable the use of shinyjs
  shinyjs::useShinyjs(),
  
  # Application title
  titlePanel("Sample size value in Network Medicine"),
  
  # Define navigation bar
  navbarPage("",
             tabPanel("Topology",
             
                      # Sidebar with a slider input for number of bins 
                      sidebarLayout(
                        sidebarPanel(
                          checkboxGroupInput("topology_dataset", label = "Dataset", 
                                             c("GTEx" = "gtex",
                                               "Scipher" = "scipher",
                                               "TCGA" = "tcga"
                                             ),
                                             selected="tcga"),
                          conditionalPanel(
                            condition = "input.topology_dataset.includes('scipher')",
                            selectInput("topology_type_dataset_scipher", label = "Type of Scipher dataset:", 
                                        choices = list("complete dataset" = "complete.dataset", "sample per patient at different times" = "sample.per.patient.all.visits", "sample per patient at baseline" = "sample.per.patient.baseline", "responders at baseline" = "responder.baseline", "non-responders at baseline" = "nonresponder.baseline"),
                                        selected = "complete.dataset"
                            ),
                          ),
                          conditionalPanel(
                            condition = "input.topology_dataset.includes('tcga')",
                            verbatimTextOutput(outputId = "topology_tcga_project"),
                            shinyWidgets::pickerInput(
                              inputId = "topology_tcga_project",
                              label = "TCGA project:",
                              choices = c("tcga", "tcga-brca", "tcga-ucec"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              selected = 'tcga',
                              multiple = TRUE),
                          ),
                          conditionalPanel(
                            condition = "input.topology_dataset.includes('gtex')",
                            verbatimTextOutput(outputId = "topology_gtex_tissues"),
                            shinyWidgets::pickerInput(
                              inputId = "topology_gtex_tissues",
                              label = "GTEx tissue:",
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
                              label = "GTEx sex:",
                              choices = c("female", "male", "both"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              selected = 'both',
                              multiple = TRUE),
                          ),
                          checkboxGroupInput("topology_method", label = "Method:", 
                                             c("ARACNE" = "aracne",
                                               "Pearson" = "pearson",
                                               "Spearman" = "spearman",
                                               "wTO" = "wto"
                                             ), selected = "pearson"
                          ),
                          checkboxGroupInput("topology_type_correlation", label = "Type of correlation:", 
                                             c("All" = "all",
                                               "Very strong" = "very strong",
                                               "Strong" = "strong",
                                               "Moderate" = "moderate",
                                               "Weak" = "weak",
                                               "Very weak" = "very weak"
                                             ),
                                             selected="all"),
                          checkboxGroupInput("topology_pvalue_threshold", label = "P-value threshold:", 
                                             choices = c("0.05" = "0.05",
                                                         "0.001" = "0.001"),
                                             selected = "0.05",
                          ),
                          selectInput("type_analysis", label = "Type analysis:",
                                      c("Topology" = "topology",
                                        "PPI" = "ppi",
                                        "Disease genes" = "disease_genes",
                                        "Essential genes" = "essential_genes"),
                                      selected = "topology"),
                          selectInput("boxplot_parameter", label = "Parameter:",
                                      c("Num. edges" = "num_edges",
                                        "Num. nodes" = "num_nodes",
                                        "Av. degree" = "av_degree",
                                        "Av. path length" = "av_path_length",
                                        "Av. clust. coef." = "av_clustering_coef",
                                        "Num. components" = "num_components",
                                        "Nodes of the LCC" = "num_lcc_nodes",
                                        "Edges of the LCC" = "num_lcc_edges",
                                        "LCC z-score" = "lcc_z",
                                        "LCC p-value" = "lcc_pvalue",
                                        "LCC log(p-value)" = "log_lcc_pvalue",
                                        "Maximum k-core" = "max_k",
                                        "Num. nodes in main core" = "num_main_core_nodes",
                                        "Num. edges in main core" = "num_main_core_edges",
                                        "Repetitions" = "rep"),
                                      selected = "num_edges"),
                          conditionalPanel(
                            condition = "input.type_analysis == 'disease_genes'",
                            verbatimTextOutput(outputId = "diseases"),
                            shinyWidgets::pickerInput(
                              inputId = "diseases",
                              label = "Diseases:",
                              choices = c("alzheimer disease", "arthritis rheumatoid", "cardiomyopathies", "diabetes mellitus type 2"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              selected = 'arthritis rheumatoid',
                              multiple = TRUE
                            )
                          ),
                          shinyWidgets::materialSwitch(inputId = "topology_log_x", label = "Turn x axis to log scale:", status="primary"),
                          shinyWidgets::materialSwitch(inputId = "topology_log_y", label = "Turn y axis to log scale:", status="primary"),
                          shinyWidgets::materialSwitch(inputId = "topology_sd", label = "Plot standard deviation:", status="primary"),
                          shinyWidgets::materialSwitch(inputId = "topology_rescale_x", label = "Normalize x axis:", status="primary"),
                          shinyWidgets::materialSwitch(inputId = "topology_rescale_y", label = "Normalize y axis:", status="primary"),
                          shinyWidgets::materialSwitch(inputId = "topology_analytical", label = "Show analytical curve:", status="primary"),
                          conditionalPanel(
                            condition = "input.topology_analytical == true",
                            #selectInput("topology_type_analytical_model", label = "Type of analytical model:", 
                            #            choices = list("exponential decay" = "exp.decay", "logarithmic" = "log", "power law" = "power.law"),
                            #            selected = "exp.decay"
                            #),
                            verbatimTextOutput(outputId = "topology_type_analytical_model"),
                            shinyWidgets::pickerInput(
                              inputId = "topology_type_analytical_model",
                              label = "Type of analytical model:",
                              choices = list("Exponential decay" = "exp.decay", "Logarithmic" = "log", "Power law" = "power.law"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              selected = 'exp.decay',
                              multiple = TRUE
                            ),
                          ),
                        ),
                        # Show boxplot of topological parameters
                        mainPanel(
                          plotOutput("topologyBoxPlot"),
                          tags$div(id = 'topologyAnalyticalModelTable')
                        )
                      )
             ),
  )
  
)