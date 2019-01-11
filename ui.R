#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

source("singleSampleAnalysisTab.R")
source("compareSampleAnalysisTab.R")
source("sidebarMenu.R")

#title <- tags$p(tags$img(src = "whitelogoi3.png", height = '30'), "DiversiTR")
#title <- tags$b("DiversiTR")

header <- dashboardHeader(titleWidth = "20%", tags$li(class = "dropdown", actionLink("resetApp", "New analysis", icon = icon("redo"))))

bodyTabs <- tabItems(
  tabItem(tabName = "aboutTab",
          fluidRow(
            box(width = 12, htmlOutput("about")),
            tags$h2("Session info"), 
            box(width = 12, verbatimTextOutput("session"))
          )),
  tabItem(tabName = "uploadRDStab",
          verbatimTextOutput("summaryRDS"),
          busyIndicator(wait = 500)
  ),
  tabItem(tabName = "uploadTXTtab",
          verbatimTextOutput("summaryTXT"),
          busyIndicator(wait = 500)
  ),
  tabItem(tabName = "showAssayTab",
          downloadButton("downloadAssay"),
          # actionButton("downloadAssay", "Download", icon = icon("download")),
          # busyIndicator(wait = 500, text = "Downloading.."),
          dataTableOutput('assayTable')
          ),
  tabItem(tabName = "showInfoTab",
          dataTableOutput("infoTable")),
  tabItem(tabName = "showMetaTab",
          #dataTableOutput("metadataTable")),
          verbatimTextOutput("metadataTable")),
  tabItem(tabName = "showHistoryTab",
          dataTableOutput("historyTable")),
  tabItem(
    tabName = "computeDiversity",
    fluidRow(
      column(
        width = 4,
        selectizeInput(
          "diversityLevel",
          "Select level",
          choices = list("V", "J", "VJ", "VpJ", "CDR3aa"),
          options = list(onInitialize = I('function() { this.setValue(""); }'))
        )
      ),
      column(
        width = 4,
        actionButton("diversityToSampleData", "Add diversity to sample data")
      )
    ),
    dataTableOutput("diversityTable"),
    htmlOutput("diversityMD")
  ),
  singleSampleTab,
  compareSampleTab
)

dashboardPage(skin = "red",
  mydashboardHeader(titleWidth = "20%", tags$li(class = "dropdown", actionLink("resetApp", "New analysis", icon = icon("refresh")))),
  dashboardSidebar(width = "20%", sideMenu),
  dashboardBody(tags$script(HTML("$('body').addClass('fixed');")),
                busyIndicator(wait = 500),
                bodyTabs
  ),
  title = "DiversiTR"
)
