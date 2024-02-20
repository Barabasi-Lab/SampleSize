library(shiny)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(gridlayout)
library(bslib)
library(dplyr)
library(DT)
library(data.table)

#shiny::runApp('/home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/SampleSizeShiny/app.R')
#shiny::runApp('/Users/j.aguirreplans/WORK/Postdoc/Projects/Scipher/SampleSize/scripts/SampleSizeShiny/app.R')


ui <- navbarPage(
  title = "MTCARS",
  selected = "Signature Visualization",
  collapsible = TRUE,
  theme = bslib::bs_theme(),
  tabPanel(
    title = "Signature Visualization",
    grid_container(
      layout = c(
        "sig_settings sig_viz"
      ),
      row_sizes = c(
        "1fr"
      ),
      col_sizes = c(
        "250px",
        "1fr"
      ),
      gap_size = "10px",
      grid_card(
        area = "sig_settings",
        card_header("Settings"),
        card_body_fill(
          div(style = "overflow-y: auto;",
              verbatimTextOutput(outputId = "input_cars"),
              shinyWidgets::pickerInput(
                inputId = "input_cars",
                label = "Select cars:",
                choices = c("Mazda RX4", "Mazda RX4 Wag", "Datsun 710"),
                options = list(
                  `actions-box` = TRUE,
                  size = 10,
                  `selected-text-format` = "count > 1",
                  container = "body"
                ),
                selected = "Mazda RX4",
                multiple = TRUE
              ),
              shinyWidgets::materialSwitch(inputId = "input_filter",
                                           label = "Filter signatures:",
                                           status = "primary"),
              conditionalPanel(
                condition = "input.input_filter == true",
                style = "margin-bottom: 20px;", # Add margin-bottom style here
                verbatimTextOutput(outputId = "input_carb"),
                shinyWidgets::pickerInput(inputId = "input_carb",
                                          label = "Filter by carb:",
                                          choices = unique(mtcars$carb),
                                          options = list(
                                            `actions-box` = TRUE,
                                            size = 10,
                                            `selected-text-format` = "count > 1"
                                          ),
                                          selected = unique(mtcars$carb),
                                          multiple = TRUE),
                verbatimTextOutput(outputId = "input_gear"),
                shinyWidgets::pickerInput(inputId = "input_gear",
                                          label = "Filter by gear:",
                                          choices = unique(mtcars$gear),
                                          options = list(
                                            `actions-box` = TRUE,
                                            size = 10,
                                            `selected-text-format` = "count > 1"
                                          ),
                                          selected = unique(mtcars$gear),
                                          multiple = TRUE),
                verbatimTextOutput(outputId = "pert_itime"),
                shinyWidgets::pickerInput(inputId = "pert_itime",
                                          label = "Filter by exposure time:",
                                          choices = list(
                                            "6 h" = "6 h",
                                            "24 h" = "24 h"
                                          ),
                                          options = list(
                                            `actions-box` = TRUE,
                                            size = 10,
                                            `selected-text-format` = "count > 1"
                                          ),
                                          multiple = TRUE),
                sliderInput(inputId = "input_mpg_cutoff",
                            label = "Filter by mpg value:",
                            value = 20,
                            min = 10,
                            max = 34),
                sliderInput(inputId = "ctl_sim_cutoff",
                            label = "Filter by control similarity:",
                            value = 0.51,
                            min = 0,
                            max = 1),
                actionButton(inputId = "filter_action",
                             label = "Filter"),
              ),
              shinyWidgets::materialSwitch(inputId = "input_merge",
                                           label = "Merge signatures:",
                                           status = "primary"),
              conditionalPanel(
                condition = "input.input_merge == true",
                verbatimTextOutput(outputId = "merge_cond"),
                shinyWidgets::pickerInput(inputId = "merge_cond",
                                          label = "Merge by same:",
                                          choices = c("carb", "gear"),
                                          options = list(
                                            `actions-box` = TRUE,
                                            size = 10,
                                            `selected-text-format` = "count > 1"
                                          ),
                                          selected = c("carb", "gear"),
                                          multiple = TRUE),
                selectInput(inputId = "merge_method",
                            label = "Method:", 
                            choices = list("DGES" = "dges",
                                           "DGES & TAS" = "dges.tas",
                                           "Spearman" = "corr"),
                            selected = "dges"
                ),
                actionButton(inputId = "merge_action",
                             label = "Merge"),
              ),
          ),
        )
      ),
      grid_card(
        area = "sig_viz",
        card_body_fill(
          tabsetPanel(
            id = "sig_viz_type",
            tabPanel(
              title = "Information",
              DTOutput(outputId = "sigTable", width = "100%")
            ),
            tabPanel(
              title = "Plot",
              grid_container(
                layout = c(
                  "sig_heatmap_area"
                ),
                row_sizes = c(
                  "1fr"
                ),
                col_sizes = c(
                  "1fr"
                ),
                gap_size = "10px",
                grid_card_plot(
                  area = "sig_heatmap_area",
                  outputId = "sigHeatmap"
                )
              )
            )
          )
        )
      )
    )
  )
)

mtcars_mod <- mtcars %>% 
  mutate(car = row.names(mtcars)) %>%
  dplyr::select(car, carb, gear, mpg)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

  #freezeReactiveValue(input, "input_cars")
  shinyWidgets::updatePickerInput(session,
                                  inputId = "input_cars",
                                  label = "Select cars:",
                                  choices = unique(mtcars_mod$car),
                                  selected = unique(mtcars_mod$car),
                                  options = list(
                                    `actions-box` = TRUE,
                                    size = 10,
                                    `selected-text-format` = "count > 1"
                                  ))
  
  process_inputs <- reactive({
    mtcars_mod_filt <- mtcars_mod %>%
      filter(car %in% input$input_cars)
    if (input$input_filter) {
      input$filter_action
      isolate({
        mtcars_mod_filt <- mtcars_mod_filt %>%
          filter((carb %in% input$input_carb) &
                   (gear %in% input$input_gear) &
                   (mpg > input$input_mpg_cutoff)
          )
      })
    }
    if (input$input_merge) {
      input$merge_action
      isolate({
        mtcars_mod_filt <- mtcars_mod_filt %>%
          group_by_at(input$merge_cond) %>%
          summarize(mpg = median(mpg)) %>%
          unite(car, all_of(input$merge_cond), sep = "|")
      })
    }
    return(mtcars_mod_filt)
  })
  
  # Observe changes in input_cars and update input_gear and input_carb choices
  observeEvent(input$input_cars, {
    selected_cars <- input$input_cars
    filtered_data <- mtcars_mod %>% filter(car %in% selected_cars)
    gear_selected <- input$input_gear
    carb_selected <- input$input_carb
    if ((is.null(gear_selected)) | (!(isTRUE(input$input_filter)))) {
      gear_selected <- unique(filtered_data$gear)
    }
    if ((is.null(carb_selected)) | (!(isTRUE(input$input_filter)))) {
      carb_selected <- unique(filtered_data$carb)
    }
    # Update input_gear choices and retain previous selections
    shinyWidgets::updatePickerInput(
      session,
      inputId = "input_gear",
      label = "Filter by gear:",
      choices = unique(filtered_data$gear),
      selected = intersect(gear_selected, unique(filtered_data$gear)),
      options = list(
        `actions-box` = TRUE,
        size = 10,
        `selected-text-format` = "count > 1"
      )
    )
    
    # Update input_carb choices and retain previous selections
    shinyWidgets::updatePickerInput(
      session,
      inputId = "input_carb",
      label = "Filter by carb:",
      choices = unique(filtered_data$carb),
      selected = intersect(carb_selected, unique(filtered_data$carb)),
      options = list(
        `actions-box` = TRUE,
        size = 10,
        `selected-text-format` = "count > 1"
      )
    )
  })
  
  output$sigTable <- DT::renderDataTable({
    mtcars_mod_filt <- process_inputs()
    mtcars_mod_filt
  })

  output$sigHeatmap <- renderPlot({
    mtcars_mod_filt <- process_inputs()
    mtcars_mod_filt %>%
      ggplot(aes(x = mpg)) +
      geom_histogram(aes(fill = car))
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
