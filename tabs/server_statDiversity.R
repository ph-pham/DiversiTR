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
        divdata <- RepSeq::divLevel(RepSeqDT(), index = meth[[1]][1], level = input$DivRepLevel, norm = normindex)
        return(data.frame(divdata, row.names = 1, check.names = FALSE))
    })
    #render plot heatmap
    output$heatmapDiv <- renderPlot({
        validate(need(!(is.null(input$DivRepLevel) || input$DivRepLevel == ""), ""))
        validate(need(!(is.null(input$DivMethod) || input$DivMethod == ""), ""))
        validate(need(!(is.null(input$DivGroup) || input$DivGroup == ""), "select group"))
        grp <- sData(RepSeqDT())[, input$DivGroup, drop=FALSE]
        #rownames(grp) <- gsub("-", ".", rownames(grp))
        dat <- as.matrix(dataDiv())
        pheatmap::pheatmap(dat, main = paste0("Diversity index: ", input$DivMethod), cluster_rows = TRUE, cluster_cols = TRUE,
            annotation_col=grp, show_colnames = T, show_rownames = T, clustering_method = "ward.D")
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