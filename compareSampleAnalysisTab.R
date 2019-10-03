source("server.R")

compareSampleTab <- tabItem(tabName = "multipleSampleTab",
                            fluidRow(
                              tabBox(width = 12,
                                id = 'multiple.topTabs',
                                tabPanel("Renyi Profile",
                                  fluidRow(
                                    # column(width = 3,
                                    #        alphaInput("")),
                                    column(width = 3,
                                      selectizeInput(
                                        "renyiLevel",
                                        "Select level",
                                        choices = list("V", "J", "VJ", "VpJ", "CDR3aa"),
                                        options = list(onInitialize = I('function() { this.setValue(""); }')))
                                    ),
                                    column(width = 3,
                                           uiOutput("renyiGroup")
                                    )#,
                                    #column(width = 6,
                                    #       uiOutput("renyiGroupChoice")
                                    #)
                                  ),  
                                  plotOutput("plotRenyi"),
                                  busyIndicator(wait = 500),
                                  htmlOutput("renyiMD"),
                                  value = "Renyitab"),

                                tabPanel(
                                  "Dissimilarity Heatmap",
                                  fluidRow(
                                    column(width = 3,
                                  selectizeInput(
                                    "dissimilarityLevel",
                                    "Select level",
                                    choices = list("V", "J", "VJ", "VpJ", "CDR3aa"),
                                    options = list(onInitialize = I('function() { this.setValue(""); }'))
                                  )),
                                  column(width = 3,
                                    selectizeInput("dissimilarityIndex", "Select dissimlarity Index",
                                        choices = list("manhattan", "euclidean", "canberra", "clark", "bray", 
                                           "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", 
                                           "mountford", "raup", "binomial", "chao", "cao", "mahalanobis"),
                                        options = list(onInitialize = I('function() { this.setValue(""); }'))
                                    )
                                  ),
                                  column(width = 3,
                                    selectizeInput("dissimilarityBinary", "Binary",
                                        choices = list("FALSE", "TRUE"),
                                        selected="FALSE") #options = list(onInitialize = I('function() { this.setValue(""); }')
                                  )
                                  ),
                                  plotOutput("plotDissimilarityHM", height = "auto", width="auto"),
                                  busyIndicator(wait = 500),
                                  withMathJax(),
                                  htmlOutput("distFuncsMD"),
                                  value = "dissimilarityHM"
                                ),
                                tabPanel(
                                  "Frequency Spectrum",
                                  fluidRow(
                                    column(width = 3,
                                           uiOutput("freqSpectrumGroup")
                                           ),
                                    column(width = 3,
                                           conditionalPanel(
                                             condition = "input.freqSpectrumGroup != null & input.freqSpectrumGroup.length>0 & input.freqSpectrumGroup != 'Sample'",
                                           radioButtons(
                                             "freqSpectrumStyle",
                                             "style du graphe",
                                             choiceNames = c("une courbe par groupe", "une courbe par echantillon"),
                                             choiceValues = c(T, F),
                                             selected = character(0)))
                                           )),
                                  plotOutput("plotfreqSpectrum"),
                                  busyIndicator(wait = 500),
                                  value = "freqSpectrum"
                                ),
                                tabPanel("Distribution clonotype (decreasing rank)",
                                         uiOutput("distribVpJGroup"),
                                         plotOutput("plotDistribVpJ"),
                                         busyIndicator(wait = 500)
                                         ),
                                tabPanel("Venn Diagram",
                                         fluidRow(
                                           column(width = 3,
                                         uiOutput("vennGroup")),
                                         column(width = 3,
                                         selectizeInput(
                                           "vennLevel",
                                           "Select level",
                                           choices = list("V", "J", "VJ", "VpJ", "CDR3aa"),
                                           options = list(onInitialize = I('function() { this.setValue(""); }'))
                                         ))),
                                         plotOutput("plotVenn", height = 600),
                                         busyIndicator(wait = 500)
                                ),
                                tabPanel("Multivariate scores",
                                         fluidRow(
                                           column(width = 3,
                                                  selectizeInput(
                                                    "muLevel",
                                                    "Select level",
                                                    choices = list("V", "J", "VJ"),
                                                    options = list(onInitialize = I('function() { this.setValue(""); }'))
                                                  )),
                                           column(width = 3,
                                                  selectizeInput(
                                                    "muType",
                                                    "Select what type",
                                                    choices = c("count", "usage"),
                                                    options = list(onInitialize = I('function() { this.setValue(""); }'))
                                                  ))
                                         ),
                                         # fluidRow(column(width = 12,
                                         #          plotOutput("plotmuScore", height = "auto"))
                                         # )
                                         plotOutput("plotmuScore", height = "auto"),
                                         busyIndicator(wait = 500)
                                  
                                ),
                                tabPanel("side by side count comparison",
                                         fluidRow(
                                           column(width = 3,
                                                  selectizeInput(
                                                    "count2v2Level",
                                                    "Select level",
                                                    choices = list("V", "J", "VJ", "VpJ", "CDR3aa"),
                                                    options = list(onInitialize = I('function() { this.setValue(""); }'))
                                                  )),
                                           column(width = 3,
                                                  uiOutput("count2v2Libs")
                                                  ),
                                           column(width = 3,
                                                  radioButtons("count2v2scale", "Choose a scale",
                                                               choices = c("counts", "log"),
                                                               selected = "counts",
                                                               inline = T)
                                                  )
                                         ),
                                         plotOutput("plot2v2count", brush = brushOpts(id = "plot2v2count_brush")),
                                         tags$hr(),
                                         uiOutput("downloadSelected"),
                                         dataTableOutput("brush2v2countDT")
                                )
                              )
                            ))
