library(shiny)
library(shinydashboard)
library(RepSeq)
library(ggplot2)
library(shinysky)
library(DT)
#source("plotFunctions.R")

options(shiny.maxRequestSize = 100 * 1024 ^ 2) # limite la taille des fichiers input, a modifier
#source("singleSampleAnalysisTab.R")

#' get data from selected area of a plot 
#' 
#' function allows to get data from selected area of a plot 
#' @param x an object of class RepSeqExperiment
#' @param level choose between V, J, VJ, VpJ, CDR3aa
#' @param libs a vector of sample names of length 2.
#' @param plot name of the previous plot
#' @return a data.table
#' @export
# @example
brush2v2count <- function(x, level = c("V", "J", "VJ", "VpJ", "CDR3aa"), libs = NULL, plot = NULL){
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    levelChoice <- match.arg(level)
    cts <- data.table::copy(assay(x))
    counts <- cts[lib ==libs[1] | lib == libs[2]]
    selected <- data.table::dcast(counts, paste(levelChoice, "~lib"), value.var="count", fun.aggregate = sum)
    bp <- shiny::brushedPoints(selected, plot)
    out <- counts[eval(parse(text = levelChoice)) %in% bp[[levelChoice]]]
    return(out)
}

alphaInput <- function(level) {
  textInput(paste0("alpha", level),
            label = "values of alpha separated by commas",
            "0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf")
}

sampleError <- function(ID){
  validate(need(!(is.null(ID) || ID == ""), "select sample"))
  
}

#Renvoie un selectizeInput
selectSample <- function(ID, sampleNames) {
  #x RepSeqExperiment object
  #id string
  selectizeInput(ID,
                 label = "Select which sample to plot",
                 choices = sampleNames,
                 options = list(onInitialize = I('function() { this.setValue(""); }')))
}
# select biological groups
selectGroup <- function(ID, x){
  sdata <- sData(x)[,unlist(lapply(sData(x), is.factor))]
  choices <- colnames(sdata)
  selectizeInput(
    ID,
    "Select group",
    choices = choices,
    options = list(onInitialize = I('function() { this.setValue(""); }'))
  )
}

#Pour plotRenyiProfile, prend str en transforme en liste de valeurs de alpha
getAlpha <- function(str) {
  as.numeric(unlist(strsplit(str, ",")))
}

convertMenuItem <- function(mi,tabName) {
  mi$children[[1]]$attribs['data-toggle']="tab"
  mi$children[[1]]$attribs['data-value'] = tabName
  mi
}


renameFiles <- function (x, y) {
  paste0(x, "/", y)
}

mydashboardHeader <- function(..., title = NULL, titleWidth = NULL, disable = FALSE,title.navbar=NULL, .list = NULL) {
  items <- c(list(...), .list)
  titleWidth <- validateCssUnit(titleWidth)
  custom_css <- NULL
  if (!is.null(titleWidth)) {
    custom_css <- tags$head(tags$style(HTML(gsub("_WIDTH_", titleWidth, fixed = TRUE, "\n      @media (min-width: 768px) {\n        .main-header > .navbar {\n          margin-left: _WIDTH_;\n        }\n        .main-header .logo {\n          width: _WIDTH_;\n        }\n      }\n    "))))
  }
  tags$header(class = "main-header", custom_css,
              style = if (disable) "display: none;",
              span(class = "logo", title),
              tags$nav(class = "navbar navbar-static-top", role = "navigation",
                       # Embed hidden icon so that we get the font-awesome dependency
                       span(shiny::icon("bars"), style = "display:none;"),
                       title.navbar,
                       div(class = "navbar-custom-menu",
                           tags$ul(class = "nav navbar-nav",
                                   items
                           )
                       )
              )
  )
}