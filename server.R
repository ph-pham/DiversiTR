#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

shinyServer(function(input, output, session) {
  # render about
  output$about<- renderUI({
    # rmarkdown::render("about.Rmd", output_format = "all")
    # knitr::knit("about.Rmd",quiet = T)
    # includeMarkdown("about.md")
    includeHTML("about.html")
    })
  # render session information
  output$session <- renderPrint({ sessionInfo() })
  # 
  output$canUpload <- reactive(return(!(
    is.null(input$chain) || input$source == ""
  )))
  # 
  outputOptions(output, "canUpload", suspendWhenHidden = FALSE)
  
  # load RDS 
  RepSeqDT <- eventReactive(c(input$samplefiles, input$RDSfile), {
    if (!is.null(input$RDSfile)) {
      RepSeqDT <- readRDS(input$RDSfile$datapath)
    } else{
      sInfo <- NULL
      if (!is.null(input$sInfofile)) { 
        sInfo <- fread(input = input$sInfofile$datapath)
        }
      tempFile <- input$samplefiles$datapath
      inFiles <- unlist(sapply(input$samplefiles$name, renameFiles, x = dirname(tempFile[1])), recursive = F)
      file.rename(tempFile, inFiles)
      RepSeqDT <- RepSeq::readClonotypeSet(
        fileList = inFiles,
        aligner = input$source,
        chain = input$chain,
        sampleinfo = sInfo
      )
      file.remove(inFiles)
    }
    return(RepSeqDT)
  })
  # output reactive
  output$isUploaded <- reactive({
    return(is.RepSeqExperiment(RepSeqDT()))
  })
  outputOptions(output, "isUploaded", suspendWhenHidden = FALSE)
  
  # open if RepSeqDT() is loaded
  observeEvent(is.RepSeqExperiment(RepSeqDT()), {
    output$showDataTab <- renderMenu({
      menuItem(
        text = "View Data in Table",
        icon = icon("table"),
        menuSubItem("show assay data",
                    tabName = "showAssayTab"),
        menuSubItem("show sample info",
                    tabName = "showInfoTab"),
        menuSubItem("show metadata",
                    tabName = "showMetaTab"),
        menuSubItem("show history",
                    tabName = "showHistoryTab"),
        menuSubItem("compute diversity indices",
                    "computeDiversity")
      )
    })
    # render single sample menu
    output$singleSampleTab <- renderMenu({
      convertMenuItem(menuItem(tabName = "singleSampleTab",
        text = "single sample analysis",
        icon = icon("user"),
        selectSample("singleSample", rownames(sData(RepSeqDT()))),
        radioButtons("singleScale", "Choose a scale",
                     choices = c("counts", "percent", "cpm"),
                     selected = character(0),
                     inline = T)
        ), "singleSampleTab")
    })
    # render multiple samples comparison menu
    output$multipleSampleTab <- renderMenu(
      menuItem(tabName = "multipleSampleTab",
        text = "comparative sample analysis",
        icon = icon("users", lib = "font-awesome")
      )
    )
    # down load RDS freshly created 
    output$downloadRDS <- renderMenu({
      menuItem("Download RDS file",
               icon = icon("download"),
               tabName = "DownloadRDS")
    })
  }) # end observeEvent
  # download RDS
  observeEvent(input$sideTabs, {
    if (input$sideTabs == "DownloadRDS") {
      dir.create(paste0(getwd(), "/", "RepSeqFiles"), showWarnings = FALSE)
      saveRDS(RepSeqDT(),
              file = paste0(getwd(), "/", "RepSeqFiles/RepSeqExperiment.rds"))
      shinyjs::info(paste0("Downloaded to ", getwd(), "/", "RepSeqFiles/RepSeqExperiment.rds"))
    }
  }, ignoreInit = T)
  
  # show summary of the RepSeqExperiment object
  output$summaryTXT <- output$summaryRDS <- renderPrint({
    #showSummary()
    flush.console()
    RepSeqDT()
  })
  # show diversity table
  showDiversity <- reactive(RepSeq::basicIndices(RepSeqDT(), level = input$diversityLevel))
  # get name
  filenameDT <- function(fname){
    return(list(list(extend = 'csv', filename = fname), list(extend = 'excel', filename = fname)))
  }
  # create download button for all assay data
  output$downloadAssay <- downloadHandler(
    "RepSeqAssay.csv",
    content = function(file) {
      write.table(assay(RepSeqDT()), file, row.names = F, sep = '\t')
    },
    contentType = "text/csv"
  )
  
  # output infoTable
  output$infoTable <- renderDataTable(sData(RepSeqDT()), 
    server = F, 
    style="bootstrap", 
    extensions = 'Buttons', 
    options = list(scrollX=TRUE, dom = 'Bfrtip', buttons = filenameDT("RepSeqInfo"))
  )
  # output AssayTable
  output$assayTable <- renderDataTable(assay(RepSeqDT()), 
    options=list(scrollX=TRUE)
  )#, server = F, style="bootstrap", extensions = 'Buttons', options = list(dom = 'Bfrtip', buttons = filenameDT("RepSeqAssay")))
  # output meta data
  #output$metadataTable <- renderDataTable(mData(RepSeqDT()), 
  #  server = F, 
  #  style="bootstrap", 
  #  extensions = 'Buttons', 
  #  options = list(scrollX=TRUE,dom = 'Bfrtip', buttons = filenameDT("RepSeqmData"))
  #)
  output$metadataTable <- renderPrint({
    out <- mData(RepSeqDT())
    if (length(out)==0) print("Nothing to display") else print(out)
  })
  # output history
  output$historyTable <- renderDataTable(History(RepSeqDT()), 
    server = F, style="bootstrap", extensions = 'Buttons', 
    options = list(scrollX=TRUE, dom = 'Bfrtip', buttons = filenameDT("RepSeqHistory"))
  )
  # output diversity 
  output$diversityTable <- renderDataTable({
    validate(need(!(is.null(input$diversityLevel) ||  input$diversityLevel == ""), "Choose a level"))
    DT::datatable(showDiversity(), style="bootstrap", extensions = 'Buttons', options = list(scrollX=TRUE, dom = 'Bfrtip', buttons = filenameDT("RepSeqDiversity")))
  })
  # plot proportion V
  output$propV <- renderPlot({
    sampleError(input$singleSample)
    validate(need(!(is.null(input$singleScale) ||  input$singleScale == ""), "Choose a scale"))
    plotPropVJ(RepSeqDT(), "V", input$singleSample, input$singleScale)
  }, height = 400)
  # plot proportion J
  output$propJ <- renderPlot({
    sampleError(input$singleSample)
    validate(need(!(is.null(input$singleScale) ||  input$singleScale == ""), "Choose a scale"))
    plotPropVJ(RepSeqDT(), "J", input$singleSample, input$singleScale)
  }, height = 400)
  
  # plot proportion VJ
  output$freqVpJ <- renderPlot({
    sampleError(input$singleSample)
    plotFreqVpJ(RepSeqDT(), input$singleSample)
  })
  # output$propJ <- output$propV <- renderPlot({
  #   if(input$tabs == "barplotJ") return(plotPropVJ(RepSeqDT(), "J", input$barplotSampleJ))
  #   if(input$tabs == "barplotV") return(plotPropVJ(RepSeqDT(), "V", input$barplotSampleV))
  # })
  # plot overlay spectratype 
  output$spectraPlot <- renderPlot({
    sampleError(input$singleSample)
    validate(need(!(is.null(input$singleScale) || input$singleScale == ""), "Choose a scale"))
    plotSpectratype(RepSeqDT(), input$singleSample, input$singleScale)
  })
  
  # plot individual spectratype 
  output$spectraPlotbis <- renderPlot({
    sampleError(input$singleSample)
    validate(need(!(is.null(input$singleScale) || input$singleScale == ""), "Choose a scale"))
    plotSpectratypebis(RepSeqDT(), input$singleSample, input$singleScale, input$spectraCDR3)
    }, width="auto", height <- function() {
    if (is.null(input$singleSample) || input$singleSample == "") return(600)
    else { 
        nrowsGrid <- ceiling(length(assay(RepSeqDT())[lib == input$singleSample, unique(V)])/4)
        return(150 * nrowsGrid)
    }
  })

  # render UI download button countheatmap
  output$downVJheatmap <- renderUI({
    if (!is.null(input$singleSample) & !(is.null(input$singleScale))) {
        downloadButton("countVJheatmap", "Download SVG")
    }
  })
  
  # plot count heatmap
  output$plotCountHeatmap <- renderPlot({
    sampleError(input$singleSample)
    validate(need(!(is.null(input$singleScale) || input$singleScale == ""), "Choose a scale"))
    plotCountVJ(x=RepSeqDT(), sampleName=input$singleSample, scale=input$singleScale)
    }, height = function() {
    if (is.null(input$singleSample) || input$singleSample == "") return(20)
    else {
      VlengthMax <- assay(RepSeqDT())[lib == input$singleSample, max(nchar(V))]
      nrowsGrid <- max(sData(RepSeqDT())$J)     
      return(30*nrowsGrid + 5*VlengthMax)
    }
  })
  
  output$countVJheatmap <- downloadHandler(
    filename =  function() {
      paste("heatmap_VJcount_", input$singleSample, ".svg", sep="")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      svg(file, height=10, width=10)
      print(plotCountVJ(x=RepSeqDT(), sampleName=input$singleSample, scale=input$singleScale)) # draw the plot
      dev.off()  # turn the device off
    },
    contentType = "image/svg"
  ) 
  
  # render plot Renyi profile
  output$renyiGroup <- renderUI({
    validate(need(!(is.null(input$renyiLevel) || input$renyiLevel ==""), "select level"))
    selectGroup("renyiGroup", RepSeqDT())
  })
  # render checkbox input
  output$renyiGroupChoice <- renderUI({
    validate(need(!(is.null(input$renyiLevel) || input$renyiLevel ==""), "select level"))
    validate(need(!(is.null(input$renyiGroup) || input$renyiGroup == ""), "select group"))
    choice <- unique(sData(RepSeqDT())[, input$renyiGroup])
    checkboxGroupInput("renyiGroupChoice", 
            label = input$renyiGroup,
            choices = choice,
            selected = choice[1],
            inline = TRUE)
  })
  # plot Renyi profile
  output$plotRenyi <- renderPlot({
    validate(need(!(is.null(input$renyiLevel) || input$renyiLevel ==""), "select level"))
    validate(need(!(is.null(input$renyiGroup) || input$renyiGroup == ""), "select group"))    
    group <- switch((input$renyiGroup == "Sample") + 1, input$renyiGroup, NULL)
    #RepSeq::plotRenyiProfiles(RepSeqDT(), c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf), input$renyiLevel, group)
    plotRenyiProfiles(x=RepSeqDT(), alpha=c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf), level=input$renyiLevel, colorBy=group)
  })  
  
  # plot dissimilarity
  output$plotDissimilarityHM <- renderPlot({
    validate(need(!(is.null(input$dissimilarityLevel) || input$dissimilarityLevel == ""), "select level"))
    validate(need(!(is.null(input$dissimilarityIndex) || input$dissimilarityIndex == ""), "select level"))
    plotDissimilarityMatrix(x=RepSeqDT(), level=input$dissimilarityLevel, method=input$dissimilarityIndex, binary=input$dissimilarityBinary)
    }, height = function() {return(as.numeric(session$clientData$output_plotDissimilarityHM_width)*0.5)}, 
        width = function() {return(as.numeric(session$clientData$output_plotDissimilarityHM_width)*0.5)})
    
  # include md
  output$distFuncsMD <- renderUI({
      shiny::withMathJax(includeMarkdown("distanceFuncs.md"))
  })
  # include md
  output$renyiMD <- renderUI({
    shiny::withMathJax(includeMarkdown("Renyi.md"))
  })
  # include md
  output$diversityMD <- renderUI({
    shiny::withMathJax(includeMarkdown("DiversityIndex.md"))
  })
  # render frequency spectrum
  output$freqSpectrumGroup <- renderUI(
    selectGroup("freqSpectrumGroup", RepSeqDT())
  )
  # plot frequency spectrum
  output$plotfreqSpectrum <- renderPlot({
    validate(need(!(is.null(input$freqSpectrumGroup) || input$freqSpectrumGroup == ""), "select group"))
    if(!input$freqSpectrumGroup == "Sample") validate(need(!(is.null(input$freqSpectrumStyle) || input$freqSpectrumStyle == ""), "select style"))
    group <- switch((input$freqSpectrumGroup == "Sample") + 1, input$freqSpectrumGroup, NULL)
    plotFrequencySpectrum(RepSeqDT(), groupBy = input$freqSpectrumStyle, colorBy = group)
  })
  
  # render VJ distribution
  output$distribVpJGroup <- renderUI(
    selectGroup("distribVpJGroup", RepSeqDT())
  )
  # plot VJ distribution
  output$plotDistribVpJ <- renderPlot({
    validate(need(!(is.null(input$distribVpJGroup) || input$distribVpJGroup == ""), "select group"))
    #group <- switch((input$distribVpJGroup == "Sample") + 1, input$distribVpJGroup, NULL)
    group <- input$distribVpJGroup
    plotDistribVpJ(RepSeqDT(), group)
  })
  # render Venn
  output$vennGroup <- renderUI(
    selectGroup("vennGroup", RepSeqDT())
  )
  # plot Venn
  output$plotVenn <- renderPlot({
    validate(need(!(is.null(input$vennLevel) || input$vennLevel == ""), "select level"))
    validate(need(!(is.null(input$vennGroup) || input$vennGroup == ""), "select group"))
    plotVenn(RepSeqDT(), level = input$vennLevel, colorBy = input$vennGroup)
  })
  
  # render mu Score 
  output$plotmuScore <- renderPlot({
    validate(need(!(is.null(input$muLevel) || input$muLevel == ""), "select level"))
    validate(need(!(is.null(input$muType) || input$muType == ""), "select type"))
    print(session$clientData)
    plotmuScore(RepSeqDT(), input$muLevel, input$muType)
    }, height = function() {    
            if(is.null(input$muLevel) || input$muLevel == "" || is.null(input$muType) || input$muType == "") return(20)
            else {
                level <- input$muLevel
                sdata <- sData(RepSeqDT())
                #groups <- sData(RepSeqDT())[,unlist(lapply(sData(RepSeqDT()), is.factor)), drop = F]
                #ngroups <- length(unique(unlist(groups[-1])))
                #nrowsGrid <- nrow(unique(assay(RepSeqDT())[,..level]))
                #return(110 + max(12*nrowsGrid, 20*ngroups + 15*ncol(groups[-1])))
                nrowsGrid <- max(sdata[, level])
                slengthMax <- max(nchar(rownames(sdata)))
                return(20 * nrowsGrid + 5 * slengthMax)
                }
        }, 
       width = function() {
            if(is.null(input$muLevel) || input$muLevel == "" || is.null(input$muType) || input$muType == "") return(20)
            else {
                level <- input$muLevel
                levelLengthMax <- assay(RepSeqDT())[, max(nchar(as.character(get(level))))]
                nsamples <- nrow(sData(RepSeqDT()))
                return(40*nsamples + 5*levelLengthMax)
                }
       }
  )
  # render 2by2 comparison
  output$count2v2Libs <- renderUI({
    selectizeInput("count2v2Libs",
                   "Select two libs",
                   choices = rownames(sData(RepSeqDT())),
                   options = list(maxItems = 2, onInitialize = I('function() { this.setValue(""); }')),
                   multiple = T
    )
  })
  # plot 2by2 comparison
  output$plot2v2count <- renderPlot({
    print(input$count2v2Libs)
    validate(need(!(is.null(input$count2v2Level) || input$count2v2Level == ""), "select level"))
    validate(need(!(is.null(input$count2v2Libs) || input$count2v2Libs == ""), "select a first sample"))
    validate(need(length(input$count2v2Libs)==2, "select a second sample"))
    plot2v2count(RepSeqDT(), input$count2v2Level, input$count2v2Libs, input$count2v2scale)
  })
  
  # selected region
  selected2v2count <- eventReactive(input$plot2v2count_brush,{
    brush2v2count(RepSeqDT(), input$count2v2Level, input$count2v2Libs, input$plot2v2count_brush)
  })
  # render selected region DT
  output$brush2v2countDT <- DT::renderDataTable({
    validate(need(!(is.null(input$count2v2Level) || input$count2v2Level == ""), "select level"))
    validate(need(!(is.null(input$count2v2Libs) || input$count2v2Libs == ""), "select a first sample"))
    validate(need(length(input$count2v2Libs)==2, "select a second sample"))
    datatable(selected2v2count(), options=list(scrollX=TRUE))
  })
  
  # dowload button for the selected data set   
  output$downloadSelected <- renderUI({
    if (!is.null(input$plot2v2count_brush)) {
        downloadButton('OutputFile', 'Download Output File')
    }
  })
  # make the download button 
  output$OutputFile <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(selected2v2count(), file, row.names=TRUE)
    }
  )
  output$singleInfoDT <- renderDataTable(sData(RepSeqDT())[input$singleSample,], options=list(scrollX=TRUE))
  output$singleAssayDT <- renderDataTable(assay(RepSeqDT())[lib == input$singleSample], options=list(scrollX=TRUE))
  
  observeEvent(input$showsampleInfo,
               showModal(modalDialog(
                 title = paste(input$singleSample," info"),
                 dataTableOutput("singleInfoDT"),
                 dataTableOutput("singleAssayDT"),
                 size = "l",
                 easyClose = T
               ))
  )
  
  observeEvent(input$resetApp, {
    session$reload()
    return()
  })
  #session$onSessionEnded(stopApp)#!! A enlever si on utilise un serveur !! stop R quand on ferme l'onglet de l'app
})