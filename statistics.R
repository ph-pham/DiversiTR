stats <- 
    tabItem("Statistical Analysis", tabName = "statisticTab",
        fluidRow(
            tabBox(width = 12,
                id = 'statA',
                tabPanel("Differential analysis",
                    fluidRow(
                        column(width = 2,
                            selectizeInput(
                                "DERepLevel",
                                "Select level",
                                choices = list("V", "J", "VJ"),
                                options = list(onInitialize = I('function() { this.setValue(""); }')))
                        ),
                        column(width = 2, 
                            uiOutput("DEGroupUI")
                        ),
                        column(width = 3,
                            uiOutput("DEcontrasts")
                        )
                    ), # end fluidRow panel Differential analysis
                    h4("PCA biplot and Volcano plot for selected choices:"), 
                    splitLayout(cellWidths = c("50%", "50%"), plotly::plotlyOutput("plotlyPCA"), plotly::plotlyOutput("volc", width = "80%")),
                    #splitLayout(cellWidths = c("50%", "50%"), plotOutput("DEplotPCA"), plotOutput("volcanoDESeq2")), 
                    hr(),
                    h4("Differential analysis results:"),
                    dataTableOutput("DETab"),
                    hr(),
                    textOutput("selectedLevel"),
                    dataTableOutput("SelectedRow"),
                    busyIndicator(wait = 500)
                ), # end panel differential analysis
                tabPanel("Diversity analysis",
                    fluidRow(
                        column(width = 2,
                            selectizeInput(
                                "DivRepLevel",
                                "Select level",
                                choices = list("V", "J", "VJ"),
                                options = list(onInitialize = I('function() { this.setValue(""); }')))
                        ),
                        column(width = 2,
                            uiOutput("DivMethodUI")
                        ),
                        column(width = 2,
                            uiOutput("DivGroupUI")
                        )
                    ), # end fluidRow panel Diversity analysis
                    plotOutput("heatmapDiv", height="auto"),
                    hr(),
                    h4("Diversity index values:"),
                    dataTableOutput("DivTab"),
                    busyIndicator(wait = 500)
                ) # end panel Diversity analysis
            ) # end tabBox
         ) # end fluidRow
    ) # end tabItems
    