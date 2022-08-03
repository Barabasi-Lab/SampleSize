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
                              choices = c("tcga", "tcga-brca", "tcga-ucec", "tcga-kirc", "tcga-lgg", "tcga-luad", "tcga-thca", "tcga-lusc", "tcga-hnsc", "tcga-prad", "tcga-skcm", "tcga-coad"),
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
                              choices = c("whole blood", "muscle - skeletal", "skin - sun exposed (lower leg)", "thyroid", "artery - tibial", "adipose - subcutaneous", "nerve - tibial", "skin - not sun exposed (suprapubic)", "lung"),
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
                                               "Very strong (>=0.8)" = "very strong",
                                               "Strong (>=0.6 & <0.8)" = "strong",
                                               "Strong or higher (>=0.6)" = "strong-very strong",
                                               "Moderate (>=0.4 & <0.6)" = "moderate",
                                               "Moderate or higher (>=0.4)" = "moderate-strong-very strong",
                                               "Weak (>=0.2 & <0.4)" = "weak",
                                               "Weak or higher (>=0.2)" = "weak-moderate-strong-very strong",
                                               "Very weak (<0.2)" = "very weak"
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
                          shinyWidgets::materialSwitch(inputId = "topology_normalize_x", label = "Normalize x axis:", status="primary"),
                          shinyWidgets::materialSwitch(inputId = "topology_normalize_y", label = "Normalize y axis:", status="primary"),
                          conditionalPanel(
                            condition = "input.topology_normalize_y == true",
                            selectInput("topology_type_normalization", label = "Type of y axis normalization:", 
                                        choices = list("Divide by max. value" = "divide.max.value", "Divide by max. possible value" = "divide.max.possible.value"),
                                        selected = "divide.max.value"
                            ),
                          ),
                          shinyWidgets::materialSwitch(inputId = "topology_analytical", label = "Show analytical curve:", status="primary"),
                          conditionalPanel(
                            condition = "input.topology_analytical == true",
                            verbatimTextOutput(outputId = "topology_type_analytical_model"),
                            shinyWidgets::pickerInput(
                              inputId = "topology_type_analytical_model",
                              label = "Type of analytical model:",
                              choices = list("Stretched exponential (by optimization)" = "Stretched exponential (by optimization)", "Stretched exponential (by linear fit)" = "Stretched exponential (by linear fit)", "Stretched exponential (without L)" = "Stretched exponential (without L)", "Logarithmic" = "Logarithmic"),
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 1"
                              ),
                              selected = 'Stretched exponential (by linear fit)',
                              multiple = TRUE
                            ),
                            selectInput("topology_type_analytical_model_output", label = "Type of output:", 
                                        choices = list("fit plot" = "fit.plot", "prediction plot" = "prediction.plot", "regression plot" = "regression.plot", "cumulative regression plot" = "cumulative.regression.plot"),
                                        selected = "fit.plot"
                            ),
                          ),
                        ),
                        # Show boxplot of topological parameters
                        mainPanel(
                          plotOutput("topologyBoxPlot", height="500px"),
                          #tableOutput('topologyAnalyticalModelTable')
                          DT::dataTableOutput('topologyAnalyticalModelTable')
                          #tags$div(id = 'topologyAnalyticalModelTable')
                          #uiOutput(outputId = "regressionPlotID", height="500px"),
                          #uiOutput(outputId = "predPlotID", height="500px"),
                        )
                      )
             ),
  )
  
)