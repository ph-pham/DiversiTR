#-------------------------------------------------------------------------------------------------------------------------------------------#
#  multiple comparison section
#-------------------------------------------------------------------------------------------------------------------------------------------#    
    # render rarefaction selection UI
    output$rareChoiceGroup <- renderUI(
        selectGroupDE("rareGroup", RepSeqDT())
    )
    # compute rarefaction table if not existing
    #raretable <- reactiveValues(raretab = NULL)
    # create reactive
    raretabreactive <- reactive({
        if (is.na(match("raretab", names(mData(RepSeqDT()))))) {
            raretab <- raretabRepSeq(RepSeqDT())
        } else {
            raretab <- mData(RepSeqDT())$raretab
        }
        return(raretab)
    })
    #observeEvent(is.RepSeqExperiment(RepSeqDT()), {
    #    raretable$raretab <- raretabreactive()
    #})
    
    # plot rarefaction curves
    output$rarecurves <- renderPlot({
        validate(need(!(is.null(input$rareGroup) || input$rareGroup == ""), ""))
        sdata <- copy(sData(RepSeqDT()))
        setDT(sdata, keep.rownames=TRUE)
        names(sdata)[1] <- "lib"
        #raretab <- raretable$raretab[sdata, on = "lib"]
        raretab <- raretabreactive()[sdata, on = "lib"]
        p <- ggplot(data = raretab, aes(x = x, y = y, fill = lib, color = get(input$rareGroup))) +
                        geom_line() + 
                        guides(fill = FALSE) + 
                        labs(title = "Rarefaction curves", 
                             x = "Number of counts", 
                             y = "Number of clonotypes", 
                             color = input$rareGroup) +
                        scale_x_continuous(breaks = pretty(1:max(raretab$x), n = 13)) + 
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                             legend.position = "top",
                             aspect.ratio = 0.5) + 
                        directlabels::geom_dl(aes(label = lib), 
                             method = list(directlabels::dl.trans(x = x + .2), "last.points"))
        return(p)
    })
    # barplot library
    output$libsizes <- plotly::renderPlotly({
    validate(need(!(is.null(input$rareGroup) || input$rareGroup == ""), ""))
        sdata <- sData(RepSeqDT())
        sdata <- data.frame(libName=rownames(sdata), sdata, check.names = FALSE)
        p <- ggplot(sdata, aes(x = libName, y = nReads)) +
            geom_bar(stat = "identity") + 
            xlab("Libraries") + ylab("Total counts") +
            theme(axis.text.x = element_text(angle = 45, hjust=1))
            #scale_x_log10(labels = scales::comma)
        p <- plotly::ggplotly(p)    
        return(p) 
    })
    
    # render downsampling choice
    output$downlibsize <- renderUI({
        sdata <- sData(RepSeqDT())
        numericInput(inputId = "libsizechoice", 
                     label = "New lib size",
                     value = min(sdata$nReads),
                     min = min(sdata$nReads),
                     max = max(sdata$nReads),
                     width = NULL
        )
    })
    
    # new libsize after downsampling # library sizes
    observeEvent(input$down, {
        output$histdownlibsizes <- renderPlot({
            cts1 <- RepSeq::assay(RepSeqDT()) 
            p1 <- histSums(cts1[, sum(count), by=VpJ][,V1], xlab="Number of counts", ylab="Number of clonotypes") + ggtitle("Clonotype Size")
            #if (is.null(RepSeqDown())) {
            #    p2 <- ggplot() +
            #            theme_void() +
            #            geom_text(aes(0,0,label='N/A')) +
            #            xlab(NULL) #optional, but safer in case another theme is applied later
            #} else {
                cts2 <- RepSeq::assay(RepSeqDown())
                p2 <- histSums(cts2[, sum(count), by=VpJ][,V1], xlab="Number of counts", ylab="Number of clonotypes") + ggtitle("Down-sampled data - Clonotype Size")
            #}
            gridExtra::grid.arrange(p1, p2, ncol=2)
        })    
    })
    # select data
    dat <- eventReactive(c(input$samplingchoice, input$down), {
        if (input$samplingchoice == "N") {
            return(RepSeqDT())
        } else {
            return(RepSeqDown())
        }
    })
    # information
    output$dataselected <- renderText({
        if (input$samplingchoice == "Y" & input$down == 1) {
            txt <- "Downsampled data set selected"
        } else {
            txt <- "Original data set selected"
        }
        return(txt)
    })
    # render plot Renyi profile
    output$renyiGroup <- renderUI({
        validate(need(!(is.null(input$renyiLevel) || input$renyiLevel ==""), "select level"))
        selectGroup("renyiGroup", dat())
    })
    # render checkbox input
    output$renyiGroupChoice <- renderUI({
        validate(need(!(is.null(input$renyiLevel) || input$renyiLevel ==""), "select level"))
        validate(need(!(is.null(input$renyiGroup) || input$renyiGroup == ""), "select group"))
        choice <- unique(sData(dat())[, input$renyiGroup])
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
        #group <- switch((input$renyiGroup == "Sample") + 1, input$renyiGroup, NULL)
        plotRenyiProfiles(x = dat(), alpha=c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf), level=input$renyiLevel, colorBy=input$renyiGroup) +
            ggplot2::theme(aspect.ratio = 0.5)
    })  
    # create Group selector
    output$GrpColMDS <- renderUI({
        validate(need(!(is.null(input$dissimilarityLevel) || input$dissimilarityLevel == ""), " "))
        validate(need(!(is.null(input$dissimilarityIndex) || input$dissimilarityIndex == ""), " "))
        selectGroup("grpCol4MDS", dat())
    })
    # plot dissimilarity
    output$plotDissimilarityHM <- renderPlot({
        validate(need(!(is.null(input$dissimilarityLevel) || input$dissimilarityLevel == ""), "select level"))
        validate(need(!(is.null(input$dissimilarityIndex) || input$dissimilarityIndex == ""), "select dissimilarity index"))
        plotDissimilarityMatrix(x = dat(), level = input$dissimilarityLevel, method = input$dissimilarityIndex, binary = FALSE)    
    })    
    # plot MDS
    output$plotMDS <- renderPlot({
        validate(need(!(is.null(input$dissimilarityLevel) || input$dissimilarityLevel == ""), " "))
        validate(need(!(is.null(input$dissimilarityIndex) || input$dissimilarityIndex == ""), " "))
        validate(need(!(is.null(input$grpCol4MDS) || input$grpCol4MDS == ""), "select group"))
        group <- switch((input$grpCol4MDS == "Sample") + 1, input$grpCol4MDS, NULL)
        plotMDS(x = dat(), level = input$dissimilarityLevel, method = input$dissimilarityIndex, colGrp = group)
    })
    # include md formula distance function
    output$distFuncsMD <- renderUI({
        shiny::withMathJax(includeMarkdown("markdown/distanceFuncs.md"))
    })
    # include md formula Reyni
    output$renyiMD <- renderUI({
        shiny::withMathJax(includeMarkdown("markdown/Renyi.md"))
    })
    # include md diversity index computation
    output$diversityMD <- renderUI({
        shiny::withMathJax(includeMarkdown("markdown/DiversityIndex.md"))
    })
    # render VJ distribution
    output$distribVpJGroup <- renderUI(
        selectGroup("distribVpJGroup", RepSeqDT())
    )
    # plot VJ distribution
    output$plotDistribVpJ <- renderPlot({
        validate(need(!(is.null(input$distribVpJGroup) || input$distribVpJGroup == ""), "select group"))
        group <- input$distribVpJGroup
        plotDistribVpJ(dat(), group, input$distribVpJGroupMeth)
    })
    # render VennUI for selecting type of Venn Diagram
    output$vennUISample <- renderUI({
    #    validate(need(!(is.null(input$vennLevel) || input$vennLevel ==""), message=" ", label=" "))
    #    validate(need(!(is.null(input$venntype) || input$venntype ==""), message=" ", label=" "))
        # render Venn samples
    #    if (input$venntype == "Samples") {
            choices <- rownames(sData(dat()))
            selectizeInput("vennSamples",
                "Select samples (3 max)",
                choices = choices,
                options = list(maxItems = 3, onInitialize = I('function() { this.setValue(""); }')),
                multiple = T)
    })
    
    #output$vennUIGroup <- renderUI({
    #    selectGroupDE("vennGroupSelected", dat())
    #})
    # render Venn
    #output$vennSubGroup <- renderUI({
    #    validate(need(!(is.null(input$vennGroupSelected) || input$vennGroupSelected == ""), "select group"))
    #    choices <- levels(sData(dat())[, input$vennGroupSelected])
    #    selectizeInput("vennSubGroup",
    #        "Select groups (3 max)",
    #        choices = choices,
    #        options = list(maxItems = 3, onInitialize = I('function() { this.setValue(""); }')),
    #        multiple = T)
    #})
    # render Venn samples
    #output$VennSamplesUI <- renderUI({
    #    validate(need(!(is.null(input$vennLevel) || input$vennLevel ==""), message=" ", label=" "))
    #    validate(need(input$vennGroupSelected == "Sample", message=" ", label=" "))
    #    choices <- rownames(sData(dat()))
    #    selectizeInput("vennSamples",
    #        "Select samples (3 max)",
    #        choices = choices,
    #        options = list(maxItems = 3, onInitialize = I('function() { this.setValue(""); }')),
    #        multiple = T
    #    )
    #})
    # plot Venn
    #output$plotVenn <- renderPlot({
    #    validate(need(!(is.null(input$vennLevel) || input$vennLevel == ""), "select level"))
    #    validate(need(!(is.null(input$venntype) || input$venntype ==""), "select type"))
    #    if( input$vennGroupSelected == "Sample") {
    #        validate(need(!(is.null(input$vennSamples) || input$vennSamples == ""), "select samples"))
    #        validate(need(length(input$vennSamples)>1, "select a second sample"))
    #        grp <- input$vennSamples 
    #    } else grp <- input$vennGroupSelected
    #        plotVenn(dat(), level = input$vennLevel, colorBy = grp)
    #})
    output$plotVenn <- renderPlot({
        validate(need(!(is.null(input$vennLevel) || input$vennLevel == ""), "select level"))
        validate(need(!(is.null(input$vennSamples) || input$vennSamples ==""), "select samples"))
        validate(need(length(input$vennSamples)>1, "select a second sample"))
        plotVenn(dat(), level = input$vennLevel, libs = input$vennSamples)
    })
    # render mu Score 
    output$plotmuScore <- renderPlot({
        validate(need(!(is.null(input$muLevel) || input$muLevel == ""), "select level"))
        validate(need(!(is.null(input$muType) || input$muType == ""), "select type"))
        print(session$clientData)
        plotmuScore(dat(), input$muLevel, input$muType)
        },  height = function() {    
                if(is.null(input$muLevel) || input$muLevel == "" || is.null(input$muType) || input$muType == "") return(20)
                else {
                    level <- input$muLevel
                    sdata <- sData(dat())
                    nrowsGrid <- max(sdata[, level])
                    slengthMax <- max(nchar(rownames(sdata)))
                    return(20 * nrowsGrid + 5 * slengthMax)
                }
            }, 
            width = function() {
                if (is.null(input$muLevel) || input$muLevel == "" || is.null(input$muType) || input$muType == "") return(20)
                else {
                    level <- input$muLevel
                    levelLengthMax <- assay(dat())[, max(nchar(as.character(get(level))))]
                    nsamples <- nrow(sData(dat()))
                    return(20*nsamples + 5*levelLengthMax)
                }
            }
    )
    # render 2by2 comparison
    output$count2v2Libs <- renderUI({
        selectizeInput("count2v2Libs",
            "Select two libs",
            choices = rownames(sData(dat())),
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
        plot2v2count(dat(), input$count2v2Level, input$count2v2Libs, input$count2v2scale) + 
            ggplot2::theme(aspect.ratio = 1)
    })
    # selected region
    selected2v2count <- eventReactive(input$plot2v2count_brush,{
        brush2v2count(dat(), input$count2v2Level, input$count2v2Libs, input$plot2v2count_brush)
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