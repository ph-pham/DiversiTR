singleSampleTab <- tabItem(tabName = "singleSampleTab",
        actionButton("showsampleInfo", "show sample Info"),
        fluidRow(
          tabBox(width = NULL,
            id = 'single.topTabs',
            tabPanel('Clonotypes frequency',
              plotOutput("freqVpJ"),
              busyIndicator(wait = 500),
              value = "freqVpJtab"
            ),
            tabPanel('Proportion',
                       plotOutput("propV"),
                       plotOutput("propJ"),
                       busyIndicator(wait = 500),
              value = "PropVJtab"
            ),
            navbarMenu('Immunoscope profile',
              tabPanel('Stacked',
                uiOutput("stackedspectraSample"),
                plotOutput("spectraPlot"),
                busyIndicator(wait = 500),
                value = "stackedspectraTypetab"
              ),
              tabPanel('Individual',
                radioButtons("spectraCDR3",
                            "Include bound CDR3 ?",
                            inline = T,
                            choiceNames = c("Yes", "No"),
                            choiceValues = c(T, F),
                            selected = F
                ),
                plotOutput("spectraPlotbis"),
                busyIndicator(wait = 500),
                value = "individualspectraTypetab"
              )
            ),
            tabPanel('VJ Distribution',
              plotOutput("plotCountHeatmap", height = "auto"),
              busyIndicator(wait = 500),
              value = "tabCountHeatmap"
            )
          )
        )
    )
