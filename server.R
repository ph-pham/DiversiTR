#--------------- shiny server ---------------#
shinyServer(function(input, output, session) {
    # render about
    output$about<- renderUI({ includeHTML("www/about.html") })
    # render session information
    output$session <- renderPrint({ sessionInfo() })
    # upload aligned & annotated CSV data 
    output$canUpload <- reactive(
        return(!(is.null(input$chain) || input$source == ""))
    )
    #
    outputOptions(output, "canUpload", suspendWhenHidden = FALSE) 
    # load RDS 
    RepSeqDT <- eventReactive(c(input$samplefiles, input$RDSfile), {
        if (!is.null(input$RDSfile)) {
        RepSeqDT <- readRDS(input$RDSfile$datapath)
        } else {
            sInfo <- NULL
            if (!is.null(input$sInfofile)) { 
                sInfo <- fread(input = input$sInfofile$datapath)
            }
            tempFile <- input$samplefiles$datapath
            inFiles <- unlist(sapply(input$samplefiles$name, renameFiles, x = dirname(tempFile[1])), recursive = F)
            file.rename(tempFile, inFiles)
            RepSeqDT <- RepSeq::readClonotypeSet(fileList = inFiles, aligner = input$source, chain = input$chain, sampleinfo = sInfo)
            file.remove(inFiles)
        }
        return(RepSeqDT)
    })
    # output reactive
    output$isUploaded <- reactive({
        return(is.RepSeqExperiment(RepSeqDT()))
    })
    outputOptions(output, "isUploaded", suspendWhenHidden = FALSE)
#-------------------------------------------------------------------------------------------------------------------------------------------#
#  Generate side menu section
#-------------------------------------------------------------------------------------------------------------------------------------------#     
source("tabs/server_sidemenu.R", local = TRUE)
    # download RDS
    observeEvent(input$sideTabs, {
        if (input$sideTabs == "DownloadRDS") {
            dir.create(paste0(getwd(), "/", "RepSeqFiles"), showWarnings = FALSE)
            saveRDS(RepSeqDT(),
            file = paste0(getwd(), "/", "RepSeqFiles/RepSeqExperiment.rds"))
            shinyjs::info(paste0("Downloaded to ", getwd(), "/", "RepSeqFiles/RepSeqExperiment.rds"))
        }
    }, ignoreInit = T)

#-------------------------------------------------------------------------------------------------------------------------------------------#
#  view RepSeqExperiment object section
#-------------------------------------------------------------------------------------------------------------------------------------------#    
source("tabs/server_viewdata.R", local = TRUE)
#-------------------------------------------------------------------------------------------------------------------------------------------#
#  single sample section
#-------------------------------------------------------------------------------------------------------------------------------------------#      
source("tabs/server_singlesample.R", local = TRUE)
#-------------------------------------------------------------------------------------------------------------------------------------------#
#  multiple comparison section
#-------------------------------------------------------------------------------------------------------------------------------------------#
source("tabs/server_multiplesamples.R", local = TRUE)
    RepSeqDown <- eventReactive(input$down, {
    #if(input$down == 0) {
    # Don't execute filter if filter-button has never been clicked.
    #    return(NULL)
    #} else {
        out <- rarefyRepSeqExp(x = RepSeqDT(), sample.size = isolate(input$libsizechoice), rngseed = isolate(input$downseed), replace = TRUE, verbose = FALSE)
        return(out)
    #}    
    })
#-------------------------------------------------------------------------------------------------------------------------------------------#
#  statististical analysis section
#-------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------#
#  diffferential analysis section
#--------------------------------------------------------------------------------#  
source("tabs/server_statDE.R", local = TRUE)
#--------------------------------------------------------------------------------#
#  diversity analysis section
#--------------------------------------------------------------------------------#
source("tabs/server_statDiversity.R", local = TRUE)
#-------------------------------------------------------------------------------------------------------------------------------------------#
#  perturbation section
#-------------------------------------------------------------------------------------------------------------------------------------------#    
source("tabs/server_perturbation.R", local = TRUE)
#-------------------------------------------------------------------------------------------------------------------------------------------#
#  download RDS section
#-------------------------------------------------------------------------------------------------------------------------------------------#    
    # sample info render
    output$singleInfoDT <- renderDataTable(sData(RepSeqDT())[input$singleSample,], options=list(scrollX=TRUE))
    # count assay render 
    output$singleAssayDT <- renderDataTable(assay(RepSeqDT())[lib == input$singleSample], options=list(scrollX=TRUE))
    # 
    observeEvent(input$showsampleInfo,
        showModal(modalDialog(
            title = paste(input$singleSample," info"),
            dataTableOutput("singleInfoDT"),
            dataTableOutput("singleAssayDT"),
            size = "l",
            easyClose = T
        ))
    )
    # reset analysis
    observeEvent(input$resetApp, {
        session$reload()
        return()
    })
  #session$onSessionEnded(stopApp)#!! comment out this commande if running on a server
}) 


