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
    output$about<- renderUI({ includeHTML("about.html") })
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
            convertMenuItem(
                menuItem(tabName = "singleSampleTab",
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
        # render statistic menu       
        output$statisticTab <- renderMenu({
            menuItem(tabName = "statisticTab",
                text = "statistical analysis",
                icon = icon("square-root-alt", lib="font-awesome")
            )
        })
        # down load RDS freshly created 
        output$downloadRDS <- renderMenu({
            menuItem("Download RDS file",
                icon = icon("download"),
                tabName = "DownloadRDS")
        })
    }) # end observeEvent RepSeqDT() is loaded
  
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
        }, contentType = "text/csv"
    )
  
    # output infoTable
    output$infoTable <- renderDataTable(RepSeq::sData(RepSeqDT()), 
        server = F, 
        style="bootstrap", 
        extensions = 'Buttons', 
        options = list(scrollX=TRUE, dom = 'Bfrtip', buttons = filenameDT("RepSeqInfo"))
    )
    # output AssayTable
    output$assayTable <- renderDataTable(RepSeq::assay(RepSeqDT()), 
        options=list(scrollX=TRUE)
    )
    # get information of slot metadata
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
    }, height = 500)
    # plot proportion J
    output$propJ <- renderPlot({
        sampleError(input$singleSample)
        validate(need(!(is.null(input$singleScale) ||  input$singleScale == ""), "Choose a scale"))
        plotPropVJ(RepSeqDT(), "J", input$singleSample, input$singleScale)
    }, height = 500)
    # plot proportion VJ
    output$freqVpJ <- renderPlot({
        sampleError(input$singleSample)
        plotFreqVpJ(RepSeqDT(), input$singleSample)
    })
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
        },  width="auto", 
            height <- function() {
                if (is.null(input$singleSample) || input$singleSample == "") return(600)
                else { 
                    nrowsGrid <- ceiling(length(assay(RepSeqDT())[lib == input$singleSample, unique(V)])/4)
                    return(150 * nrowsGrid)
                }
            }
    )
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
        },  height = function() {
                if (is.null(input$singleSample) || input$singleSample == "") return(20)
                else {
                    VlengthMax <- assay(RepSeqDT())[lib == input$singleSample, max(nchar(V))]
                    nrowsGrid <- max(sData(RepSeqDT())$J)     
                    return(30*nrowsGrid + 5*VlengthMax)
                }
            }
    )
    # download button for count heatmap
    output$countVJheatmap <- downloadHandler(
        filename =  function() {
            paste("heatmap_VJcount_", input$singleSample, ".svg", sep="")
        },
        # content is a function with argument file. content writes the plot to the device
        content = function(file) {
            svg(file, height=10, width=10)
            print(plotCountVJ(x=RepSeqDT(), sampleName=input$singleSample, scale=input$singleScale)) # draw the plot
            dev.off()  # turn the device off
        }, contentType = "image/svg"
    )
#-------------------------------------------------------------------------------------------------------------------------------------------#
#  multiple comparison section
#-------------------------------------------------------------------------------------------------------------------------------------------#    
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
        validate(need(!(is.null(input$renyiLevel) || input$renyiLevel ==""), " "))
        validate(need(!(is.null(input$renyiGroup) || input$renyiGroup == ""), " "))    
        group <- switch((input$renyiGroup == "Sample") + 1, input$renyiGroup, NULL)
        plotRenyiProfiles(x=RepSeqDT(), alpha=c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf), level=input$renyiLevel, colorBy=group)
    })  
    # plot dissimilarity
    output$plotDissimilarityHM <- renderPlot({
        validate(need(!(is.null(input$dissimilarityLevel) || input$dissimilarityLevel == ""), "select level"))
        validate(need(!(is.null(input$dissimilarityIndex) || input$dissimilarityIndex == ""), "select dissimilarity index"))
        plotDissimilarityMatrix(x=RepSeqDT(), level=input$dissimilarityLevel, method=input$dissimilarityIndex, binary="FALSE")    
    })
    # create Group selector
    output$GrpColMDS <- renderUI({
        validate(need(!(is.null(input$dissimilarityLevel) || input$dissimilarityLevel == ""), " "))
        validate(need(!(is.null(input$dissimilarityIndex) || input$dissimilarityIndex == ""), " "))
        selectGroup("grpCol4MDS", RepSeqDT())
    })
    # plot MDS  
    output$plotMDS <- renderPlot({
        validate(need(!(is.null(input$dissimilarityLevel) || input$dissimilarityLevel == ""), " "))
        validate(need(!(is.null(input$dissimilarityIndex) || input$dissimilarityIndex == ""), " "))
        validate(need(!(is.null(input$grpCol4MDS) || input$grpCol4MDS == ""), "select group"))
        group <- switch((input$grpCol4MDS == "Sample") + 1, input$grpCol4MDS, NULL)
        plotMDS(RepSeqDT(), level=input$dissimilarityLevel, method=input$dissimilarityIndex, colGrp=group)
    })
    # include md formula distance function
    output$distFuncsMD <- renderUI({
        shiny::withMathJax(includeMarkdown("distanceFuncs.md"))
    })
    # include md formula Reyni
    output$renyiMD <- renderUI({
        shiny::withMathJax(includeMarkdown("Renyi.md"))
    })
    # include md diversity index computation
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
        group <- input$distribVpJGroup
        plotDistribVpJ(RepSeqDT(), group, input$distribVpJGroupMeth)
    })
    # render Venn
    output$vennGroup <- renderUI(
        selectGroup("vennGroupSelected", RepSeqDT())
    )
    # render Venn samples
    output$VennSamplesUI <- renderUI({
        validate(need(!(is.null(input$vennLevel) || input$vennLevel ==""), message=" ", label=" "))
        validate(need(input$vennGroupSelected == "Sample", message=" ", label=" "))
        choices <- rownames(sData(RepSeqDT()))
        selectizeInput("vennSamples",
            "Select samples (3 max)",
            choices = choices,
            options = list(maxItems = 3, onInitialize = I('function() { this.setValue(""); }')),
            multiple = T
        )
    })
    # plot Venn
    output$plotVenn <- renderPlot({
        validate(need(!(is.null(input$vennLevel) || input$vennLevel == ""), "select level"))
        validate(need(!(is.null(input$vennGroupSelected) || input$vennGroupSelected == ""), "select group"))
        if( input$vennGroupSelected == "Sample") {
            validate(need(!(is.null(input$vennSamples) || input$vennSamples == ""), "select samples"))
            validate(need(length(input$vennSamples)>1, "select a second sample"))
            grp <- input$vennSamples 
        } else grp <- input$vennGroupSelected
            plotVenn(RepSeqDT(), level = input$vennLevel, colorBy = grp)
    })
    # render mu Score 
    output$plotmuScore <- renderPlot({
        validate(need(!(is.null(input$muLevel) || input$muLevel == ""), "select level"))
        validate(need(!(is.null(input$muType) || input$muType == ""), "select type"))
        print(session$clientData)
        plotmuScore(RepSeqDT(), input$muLevel, input$muType)
        },  height = function() {    
                if(is.null(input$muLevel) || input$muLevel == "" || is.null(input$muType) || input$muType == "") return(20)
                else {
                    level <- input$muLevel
                    sdata <- sData(RepSeqDT())
                    nrowsGrid <- max(sdata[, level])
                    slengthMax <- max(nchar(rownames(sdata)))
                    return(20 * nrowsGrid + 5 * slengthMax)
                }
            }, 
            width = function() {
                if (is.null(input$muLevel) || input$muLevel == "" || is.null(input$muType) || input$muType == "") return(20)
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
#-------------------------------------------------------------------------------------------------------------------------------------------#
#  statististical analysis section
#-------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------#
#  diffferential analysis section
#--------------------------------------------------------------------------------#  
    # render select group UI
    output$DEGroupUI <- renderUI({
        validate(need(!(is.null(input$DERepLevel) || input$DERepLevel == ""), "select level"))
        selectGroup("DEGroup", RepSeqDT())
    })
    # create data DESEq2 
    dataDESeq <- reactive({
        validate(need(!(is.null(input$DERepLevel) || input$DERepLevel == ""), ""))
        validate(need(!(is.null(input$DEGroup) || input$DEGroup == ""), ""))
        datatabDESeq2 <- RepSeq::toDESeq2(RepSeqDT(), conditions=input$DEGroup, level=input$DERepLevel)
        datatabDESeq2 <- DESeq2::estimateSizeFactors(datatabDESeq2, type="poscounts")
        datatabDESeq2 <- DESeq2::DESeq(datatabDESeq2, fitType='local')
        return(datatabDESeq2)
    })
    # plot normalized PCA (not run)
    output$DEplotPCA <- renderPlot({      
        vsd <- DESeq2::rlog(dataDESeq())
        datapca <- DESeq2::plotPCA(vsd, intgroup=input$DEGroup, returnData = TRUE)
        percentVar <- round(100 * attr(datapca, "percentVar"))
        p <- ggplot2::ggplot(datapca, ggplot2::aes_string(x="PC1", y="PC2", color=input$DEGroup)) + ggplot2::geom_point(size=4) + 
                ggplot2::xlab(paste0("PC1: ", percentVar[1], "% variance")) +
                ggplot2::ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
                ggplot2::ggtitle(paste0("PCA Biplot for ", input$DEGroup)) +
                ggplot2::coord_fixed(ratio=2) + 
                ggplot2::theme(axis.text = ggplot2::element_text(size = 16), axis.title = ggplot2::element_text(size = 18))
        print(p)
    })
    output$plotlyPCA <- plotly::renderPlotly({
        vsd <- DESeq2::rlog(dataDESeq())
        datapca <- DESeq2::plotPCA(vsd, intgroup=input$DEGroup, returnData = TRUE)
        percentVar <- round(100 * attr(datapca, "percentVar"))
        p <- ggplot2::ggplot(datapca, ggplot2::aes_string(x="PC1", y="PC2", color=input$DEGroup)) + ggplot2::geom_point(aes(text=rownames(datapca)), size=4) + 
                ggplot2::xlab(paste0("PC1: ", percentVar[1], "% variance")) +
                ggplot2::ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
                ggplot2::ggtitle(paste0("PCA Biplot for ", input$DEGroup)) +
                ggplot2::coord_fixed(ratio=2) + 
                ggplot2::theme(axis.text = ggplot2::element_text(size = 16), axis.title = ggplot2::element_text(size = 18))
        p <- plotly::ggplotly(p)
    })
    # render 2by2 comparison
    output$DEcontrasts <- renderUI({
        validate(need(!(is.null(input$DERepLevel) || input$DERepLevel == ""), ""))
        validate(need(!(is.null(input$DEGroup) || input$DEGroup == ""), ""))
        selectizeInput("GrpSelect",
            "Comparison",
            choices = levels(SummarizedExperiment::colData(dataDESeq())[, input$DEGroup]),
            options = list(maxItems = 2, onInitialize = I('function() { this.setValue(""); }')),
            multiple = T
        )
    })
    # create DESeq2 result tab
    DESeqRes <- reactive({
        validate(need(!(is.null(input$DERepLevel) || input$DERepLevel == ""), ""))
        validate(need(!(is.null(input$DEGroup) || input$DEGroup == ""), ""))
        validate(need(!(is.null(input$GrpSelect) || input$GrpSelect == ""), "select first group"))
        validate(need(length(input$GrpSelect)==2, "select second group"))
        res <- DESeq2::results(dataDESeq(), contrast=c(input$DEGroup, input$GrpSelect))
        res <- res[order(res$padj),]    
        return(as.data.frame(res))
    })
    # create volcano plot not used
    output$volcanoDESeq2 <- renderPlot({
        tmp <- DESeqRes()
        tmp <- tmp[!is.na(tmp$padj), ]
        with(tmp, plot(log2FoldChange, -log10(padj), pch=20, main=paste0("Volcano plot: ", input$DEGroup), cex=1.0, 
                xlab=bquote(~Log[2]~fold~change), 
                ylab=bquote(~-log[10]~adj~pvalue)))
        with(subset(tmp, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
        #Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
        abline(v=0, col="black", lty=3, lwd=1.0)
        abline(v=-2, col="black", lty=4, lwd=2.0)
        abline(v=2, col="black", lty=4, lwd=2.0)
        if (length(tmp$pvalue[tmp$padj<0.05])>0) {
            abline(h=-log10(max(tmp$pvalue[tmp$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
            gn.selected <- abs(tmp$log2FoldChange) > 2.5 & tmp$padj < 0.01
            if (any(gn.selected)) text(tmp$log2FoldChange[gn.selected], -log10(tmp$padj)[gn.selected], lab=rownames(tmp)[gn.selected], cex=0.6, pos=2)
        }
    })
    # volcano plot plotly
    output$volc <- plotly::renderPlotly({
        degTab <- DESeqRes()
        degTab <- degTab[!is.na(degTab$padj), ]
        data.table::setDT(degTab, keep.rownames=TRUE)
        degTab[, group := "NotSignificant"]
        # change the grouping for the entries with both significance and large enough fold change
        degTab[padj < 0.05 & abs(log2FoldChange) >= 2, group := "Significant&FoldChange"]
        # change the grouping for the entries with significance but not a large enough Fold change
        degTab[padj < 0.05 & abs(log2FoldChange) < 2, group := "Significant"]
        # change the grouping for the entries a large enough Fold change but not a low enough p value
        degTab[padj > 0.05 & abs(log2FoldChange) >= 2, group := "FoldChange"]
        degTab[, BHpvalue := -log10(padj)]
        p <- degTab[, plotly::plot_ly(x = log2FoldChange, y = BHpvalue, text = rn, mode = "markers", color = group)] 
        p <- p %>% plotly::layout(title ="Volcano Plot") %>% plotly::layout(legend = list(orientation = 'h')) %>% plotly::toWebGL()
        p$elementId <- NULL
        print(p)
    })
    # render DESeq2 results 
    output$DETab <- renderDataTable({
        return(datatable(DESeqRes(), options = list(scrollX=TRUE, dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel'), pageLength = 10)) %>% formatRound(c(1:4), 2) %>% formatRound(c(5:6), 4))
        }, server=FALSE
    )
    # render text 
    output$selectedLevel <- renderText({ 
        validate(need(!(is.null(input$GrpSelect) || input$GrpSelect == ""), ""))
        validate(need(length(input$GrpSelect)==2, ""))    
        s <- rownames(DESeqRes())[input$DETab_rows_selected]
        paste("Assay data for combination of level", input$DERepLevel, ":", s)
    })
    # get data for selected row
    output$SelectedRow <- renderDataTable({
        validate(need(!(is.null(input$GrpSelect) || input$GrpSelect == ""), ""))
        validate(need(length(input$GrpSelect)==2, ""))
        validate(need(!(is.null(input$DETab_rows_selected) || input$DETab_rows_selected == ""), "select row"))
        s <- rownames(DESeqRes())[input$DETab_rows_selected]
        dcast(RepSeq::assay(RepSeqDT())[get(input$DERepLevel) %in% s], VpJ~lib, value.var="count", fun.aggregate=sum)
        }, style="bootstrap", extensions = 'Buttons', options=list(scrollX=TRUE, dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel'))
    )
#--------------------------------------------------------------------------------#
#  diversity analysis section
#--------------------------------------------------------------------------------#    
    # render UI choice of diversity index
    output$DivMethodUI <- renderUI({
        validate(need(!(is.null(input$DivRepLevel) || input$DivRepLevel == ""), "select level"))
        selectizeInput("DivMethod",
            "Select method",
            choices = list("shannon", "shannon.norm", "simpson", "simpson.norm", "invsimpson", "inv.norm") #,
            #options = list(onInitialize = I('function() { this.setValue(""); }'))
        )
    })
    # render select gropu UI for diversity analysis
    output$DivGroupUI <- renderUI({
        validate(need(!(is.null(input$DivRepLevel) || input$DivRepLevel == ""), ""))
        validate(need(!(is.null(input$DivMethod) || input$DivMethod == ""), "select method"))
        selectGroup("DivGroup", RepSeqDT())
    })
    # compute diversity at a level of th repertoire
    dataDiv <- reactive({
        validate(need(!(is.null(input$DivRepLevel) || input$DivRepLevel == ""), ""))
        validate(need(!(is.null(input$DivMethod) || input$DivMethod == ""), ""))
        meth <- strsplit(input$DivMethod, split="\\.")
        if (is.na(meth[[1]][2])) normindex <- TRUE else normindex <- FALSE
        divdata <- RepSeq::divLevel(RepSeqDT(), index=meth[[1]][1], level=input$DivRepLevel,norm=normindex)
        return(data.frame(divdata, row.names=1))
    })
    #render plot heatmap
    output$heatmapDiv <- renderPlot({
        validate(need(!(is.null(input$DivRepLevel) || input$DivRepLevel == ""), ""))
        validate(need(!(is.null(input$DivMethod) || input$DivMethod == ""), ""))
        validate(need(!(is.null(input$DivGroup) || input$DivGroup == ""), "select group"))
        grp <- sData(RepSeqDT())[, input$DivGroup, drop=FALSE]
        rownames(grp) <- gsub("-", ".", rownames(grp))
        dat <- as.matrix(dataDiv())
        pheatmap::pheatmap(dat, main=paste0("Diversity index: ", input$DivMethod), cluster_rows = TRUE, cluster_cols = TRUE,
            annotation_col=grp, show_colnames=T, show_rownames=T, clustering_method = "ward.D")
        }, height = function() {    
                nrowsGrid <- nrow(dataDiv())
                slengthMax <- max(nchar(rownames(dataDiv())))
                return(20 * nrowsGrid + 5 * slengthMax)
           }
    )
    # render table of diversity index computed at le chosen level
    output$DivTab <- renderDataTable({
        return(datatable(dataDiv(), options = list(scrollX=TRUE, dom = 'Bfrtip', pageLength = 10)) %>% formatRound(c(1: ncol(dataDiv())), 4))
    })
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
  #session$onSessionEnded(stopApp)#!! A enlever si on utilise un serveur !! stop R quand on ferme l'onglet de l'app
}) 


