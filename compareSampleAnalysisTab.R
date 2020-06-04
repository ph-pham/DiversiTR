#source("server.R")

compareSampleTab <- 
    tabItem(tabName = "multipleSampleTab",
        fluidRow(
            tabBox(width = 12,
                id = 'multiple.topTabs',
                tabPanel("Renyi Profile",
                    fluidRow(
                        column(width = 3,
                            selectizeInput(
                                "renyiLevel",
                                "Select level",
                                choices = list("V", "J", "VJ", "VpJ", "CDR3aa"),
                                options = list(onInitialize = I('function() { this.setValue(""); }')))
                        ),
                        column(width = 3,
                            uiOutput("renyiGroup")
                        )
                    ),  
                    plotOutput("plotRenyi"),
                    busyIndicator(wait = 500),
                    htmlOutput("renyiMD"),
                    value = "Renyitab"
                ),
                tabPanel("Dissimilarity Analysis",
                    fluidRow(
                        column(width = 3,
                            selectizeInput(
                                "dissimilarityLevel",
                                "Select level",
                                choices = list("V", "J", "VJ", "VpJ", "CDR3aa"),
                                options = list(onInitialize = I('function() { this.setValue(""); }'))
                            )
                        ),
                        column(width = 3,
                            selectizeInput("dissimilarityIndex", "Select dissimlarity Index",
                                choices = list("manhattan", "euclidean", "canberra", "clark", "bray", 
                                    "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", 
                                    "mountford", "raup", "binomial", "chao", "cao", "mahalanobis"),
                                options = list(onInitialize = I('function() { this.setValue(""); }'))
                            )
                        ),
                        column(width = 3,
                            uiOutput("GrpColMDS")
                        )
                    ),
                    splitLayout(cellWidths = c("50%", "50%"), plotOutput("plotDissimilarityHM"), plotOutput("plotMDS")),
                        busyIndicator(wait = 500),
                        withMathJax(),
                        htmlOutput("distFuncsMD"),
                        value = "dissimilarityHM"
                ),
                tabPanel("Frequency Spectrum",
                    fluidRow(
                        column(width = 3,
                            uiOutput("freqSpectrumGroup")
                        ),
                        column(width = 3,
                            conditionalPanel(
                                condition = "input.freqSpectrumGroup != null & input.freqSpectrumGroup.length>0 & input.freqSpectrumGroup != 'Sample'",
                                radioButtons("freqSpectrumStyle",
                                    "Graph format",
                                    choiceNames = c("Curv by group", "Curv by sample"),
                                    choiceValues = c(T, F),
                                    selected = character(0)
                                )
                            )
                        )
                    ),
                    plotOutput("plotfreqSpectrum"),
                    busyIndicator(wait = 500),
                    value = "freqSpectrum"
                ),
                tabPanel("Distribution clonotype",
                    fluidRow(
                        column(width = 3, 
                            uiOutput("distribVpJGroup")
                        ),
                        column(width = 3, 
                            selectInput("distribVpJGroupMeth", 
                                "Select method",
                                choices = c("Sum" = "sum", "Average" = "mean"),
                                selected = "Sum"
                            )
                        )
                    ),
                    plotOutput("plotDistribVpJ"),
                    busyIndicator(wait = 500)
                ),
                tabPanel("Venn Diagram",
                    fluidRow(
                        column(width = 3,
                            uiOutput("vennGroup")
                        ),
                        column(width = 3,
                            selectizeInput("vennLevel",
                                "Select level",
                                choices = list("V", "J", "VJ", "VpJ", "CDR3aa"),
                                options = list(onInitialize = I('function() { this.setValue(""); }'))
                            )
                        ),
                        column(width = 3,
                            uiOutput("VennSamplesUI")
                        )
                    ),
                    plotOutput("plotVenn", height = 800, width = 800),
                    busyIndicator(wait = 500)
                ),
                tabPanel("Multivariate scores",
                    fluidRow(
                        column(width = 3,
                            selectizeInput("muLevel",
                                "Select level",
                                choices = list("V", "J", "VJ"),
                                options = list(onInitialize = I('function() { this.setValue(""); }'))
                            )
                        ),
                        column(width = 3,
                            selectizeInput("muType",
                                "Select data type",
                                choices = c("count", "usage"),
                                options = list(onInitialize = I('function() { this.setValue(""); }'))
                            )
                        )
                    ),
                    plotOutput("plotmuScore", height = "auto"),
                    busyIndicator(wait = 500)                            
                ),
                tabPanel("Two samples comparison",
                    fluidRow(
                        column(width = 3,
                            selectizeInput("count2v2Level",
                                "Select level",
                                choices = list("V", "J", "VJ", "VpJ", "CDR3aa"),
                                options = list(onInitialize = I('function() { this.setValue(""); }'))
                            )
                        ),
                        column(width = 3,
                            uiOutput("count2v2Libs")
                        ),
                        column(width = 3,
                            radioButtons("count2v2scale", 
                                "Choose a scale",
                                choices = c("counts", "log"),
                                selected = "counts",
                                inline = T
                            )
                        )
                    ),
                    plotOutput("plot2v2count", brush = brushOpts(id = "plot2v2count_brush")),
                    tags$hr(),
                    uiOutput("downloadSelected"),
                    dataTableOutput("brush2v2countDT")
                )
            )
        )
    )
