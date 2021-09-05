require(waiter)
require(DT)
require(ggplot2)
require(cowplot)
require(GA)

library(shiny)

waiting_screen <- tagList(
    spin_solar(),
    HTML("<br><br><br><br><br>Generating model output...")
)

ui <- bootstrapPage(
    mainPanel(
        # ui <- fluidPage(
        use_waiter(),
        # Application title
        titlePanel("U-series isotopes gain loss model - genetic algorithm optimisation"),

        tabsetPanel(id = "abset",
                    # Load data panel
                    tabPanel("Load the data",
                             p("Before uploading, check that your CSV file contains columns with these names:"),
                             HTML("
                            <li> <b>Sample ID</b>: sample identification code
                            <li> <b>Depth (m)</b>: sample depth (in m)
                            <li> <b>RU234_U238</b>: [<sup>234</sup>U/<sup>238</sup>U] activity ratios
                            <li> <b>R48_2SE</b>: the 2 sigma errors of the [<sup>234</sup>U/<sup>238</sup>U] activity ratios
                            <li> <b>RTh230_U238</b>: [<sup>230</sup>Th/<sup>238</sup>U] activity ratios
                            <li> <b>R08_2SE</b>: the 2 sigma errors of the [<sup>230</sup>Th/<sup>238</sup>U] activity ratios <br><br>
                            The deepest sample is taken at the starting composition, so make sure that what you want to use as starting material is included and assigned the deepest depth <br><br>
                            Please cite:<br> <b>Dosseto, A., et al. (2014). 'Age and weathering rate of sediments in small catchments: The role of hillslope erosion.' Geochimica et Cosmochimica Acta 132(0): 238-258 </b>
                            "),
                             tags$hr(),
                             fileInput("file1",
                                       "Choose CSV file",
                                       accept = c("text/csv",
                                                  "text/comma-separated-values,text/plain",
                                                  ".csv")),
                             actionButton('gotoinspect', 'Go to inspect the data')
                    ), # end of tab

                    # Inspect data panel
                    tabPanel("Inspect the data",
                             value = "inspectthedata",
                             p("Here is the raw data from the CSV file"),
                             DT::dataTableOutput('contents'),
                             actionButton('gotosetmodel', 'Go to set the model parameters')
                    ), # end of tab

                    # Set model parameters panel
                    tabPanel("Set model parameters",
                             value = "setmodelparameters",
                             fluidRow(
                                 column(4,
                                        numericInput("nbit", "Number of iterations:", 1000, min = 1, max = 1e6),
                                        numericInput("logT_min", "Minimum log10 of weathering age of top sample (in log(yr)):", 2, min = 1, max = 7),
                                        numericInput("logT_max", "Maximum log10 of weathering age of top sample (in log(yr)):", 6, min = 1, max = 7),
                                        numericInput("Tmult_min", "Minimum fraction weathering age of sample compared to that of sample directly above:", 0.1, min = 0.001, max = 1),
                                 ),
                                 column(4,
                                        numericInput("logk238_min", "Minimum log10 of k238 (in log(yr-1)):", -8, min = -20, max = 0),
                                        numericInput("logk238_max", "Maximum log10 of k238 (in log(yr-1)):", -3, min = -20, max = 0),
                                        numericInput("k48_min", "Minimum k234/k238:", 0.1, min = 0.01, max = 1000),
                                        numericInput("k48_max", "Maximum k234/k238:", 10, min = 0.01, max = 1000)
                                 ),
                                 column(4,
                                        numericInput("logf238_k238_min", "Minimum log10 of f238/k238:", -8, min = -20, max = 3),
                                        numericInput("logf238_k238_max", "Maximum log10 of f238/k238:", 3, min = -20, max = 3),
                                        numericInput("f48_min", "Minimum f234/f238:", 0.1, min = 0.1, max = 1000),
                                        numericInput("f48_max", "Maximum f234/f238:", 10, min = 0.1, max = 1000),
                                        # save solutions
                                        textInput(inputId = "filename", label = "File name to save calculated ratios and weathering ages to", value = "profile_name")
                                 )),
                             actionButton("run", label = "Run simulation and visualise the output")
                    ),

                    # Visualise results tab
                    tabPanel("Visualise the model",
                             value = "visualise",
                             tags$hr(),
                             textOutput("print_results"),
                             tags$hr(),
                             textOutput("param_sol"),
                             tableOutput("print_param"),
                             # show a spinner while we wait for the plots to draw
                             #  withSpinner()
                             # # Input: Choose dataset ----
                             # selectInput("dataset", "Choose a dataset:",
                             #             choices = c("rock", "pressure", "cars")),
                             #
                             # Button
                             downloadButton("downloadData", "Download"),
                             tags$hr(),
                             plotOutput("plots",
                                        width = "100%",
                                        height = "300px"),
                             color="blue",
                             size = 5
                    ) # end of tab

        )
    ))

# Define server logic
server <- function(input, output, session) {
    # activate the buttons to move between tabs
    observeEvent(input$gotoinspect, {
        updateTabsetPanel(session, "abset",
                          selected = "inspectthedata")
    })

    output$contents <- DT::renderDataTable({
        inFile <- input$file1

        if (is.null(inFile))
            return(DF)

        read.csv(inFile$datapath)
    })

    observeEvent(input$gotosetmodel, {
        updateTabsetPanel(session, "abset",
                          selected = "setmodelparameters"  )
    })

    observeEvent(input$run, {
        updateTabsetPanel(session, "abset",
                          selected = "visualise"  )
    })


    observeEvent(input$run, {
        waiter_show(
            html = waiting_screen,
            color = "#333e48",
            id = 'plots',
            hide_on_render = TRUE
        )

        showNotification("Starting model run...")

        # output_model <-  reactive({
        inFile <- input$file1
        if (is.null(inFile)) {
            input_data <- DF
        } else {
            input_data <- read.csv(inFile$datapath)
        }

        output_model <- UseriesGAwigl::U_series_gain_loss_function(input_data, nbit = input$nbit,
                                                    logT_min = input$logT_min, logT_max = input$logT_max,
                                                    Tmult_min = input$Tmult_min,
                                                    logk238_min = input$logk238_min, logk238_max = input$logk238_max,
                                                    k48_min = input$k48_min, k48_max = input$k48_max,
                                                    logf238_k238_min = input$logf238_k238_min, logf238_k238_max = input$logf238_k238_max,
                                                    f48_min = input$f48_min, f48_max = input$f48_max)
        # output_model
        # })

        showNotification("Model run complete.")


        output$print_results <- renderText({ paste("Regolith production rate (mm/kyr):",
                                                   round(output_model$regolith_production_rate,1)) })

        output$param_sol <- renderText("Solutions:")

        output_model$parameters[c(1,3)] <- log10(output_model$parameters[c(1,3)])


        names(output_model$parameters) <- c("log10(k238)", "k234/k238", "log10(f238/k238)",
                                            "f234/f238", "r48i", "r08i")

        output$print_param <- renderTable({ output_model$parameters },
                                          colnames = TRUE)

        names(output_model$calculated_ratios_and_ages) <- c("[234U/238U]", "[230Th/238U]",
                                                            "Weathering age (yr)", "Depth (m)")

        # save solutions for weathering ages
        output$downloadData <- downloadHandler(
            filename = function() {
                paste(input$filename, ".csv", sep = "")
            },
            content = function(file) {
                write.csv(output_model$calculated_ratios_and_ages, file, row.names = F)
            }
        )

        # # draw the plots
        output$plots <- renderPlot({
            # print(output_model$graph_t) })
            #     # draw plots in a panel
            print(output_model$graphs)
        })
        # })
    })
}

# Run the application
shinyApp(ui = ui, server = server)
