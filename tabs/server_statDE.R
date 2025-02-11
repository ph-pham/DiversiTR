#--------------------------------------------------------------------------------#
#  statistical analysis - diffferential analysis section
#--------------------------------------------------------------------------------#  
    # render select group UI
    output$DEGroupUI <- renderUI({
        validate(need(!(is.null(input$DERepLevel) || input$DERepLevel == ""), "select level"))
        selectGroupDE("DEGroup", RepSeqDT())
    })
    # create data DESEq2 
    dataDESeq <- reactive({
        validate(need(!(is.null(input$DERepLevel) || input$DERepLevel == ""), ""))
        validate(need(!(is.null(input$DEGroup) || input$DEGroup == ""), ""))
        datatabDESeq2 <- RepSeq::toDESeq2(RepSeqDT(), conditions=input$DEGroup, level=input$DERepLevel)
        datatabDESeq2 <- DESeq2::estimateSizeFactors(datatabDESeq2, type="poscounts")
        datatabDESeq2 <- DESeq2::DESeq(datatabDESeq2, fitType="local")
        return(datatabDESeq2)
    })
    # plot normalized PCA (not run)
    output$DEplotPCA <- renderPlot({
        #validate(need(!(is.null(input$DERepLevel) || input$DERepLevel == ""), ""))
        #validate(need(!(is.null(input$DEGroup) || input$DEGroup == ""), ""))    
        vsd <- DESeq2::varianceStabilizingTransformation(dataDESeq())
        datapca <- DESeq2::plotPCA(vsd, intgroup=input$DEGroup, returnData = TRUE)
        percentVar <- round(100 * attr(datapca, "percentVar"))
        p <- ggplot2::ggplot(datapca, ggplot2::aes_string(x="PC1", y="PC2", color=input$DEGroup)) + ggplot2::geom_point(size=4) + 
                ggplot2::xlab(paste0("PC1: ", percentVar[1], "% variance")) +
                ggplot2::ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
                ggplot2::ggtitle("Principal Component Analysis") +
                ggplot2::coord_fixed(ratio=2) + 
                ggplot2::theme(axis.text = ggplot2::element_text(size = 16), axis.title = ggplot2::element_text(size = 18))
        print(p)
    })
    # plotly PCA performed on normalised data (method vst).
    output$plotlyPCA <- plotly::renderPlotly({
        validate(need(!(is.null(input$DERepLevel) || input$DERepLevel == ""), ""))
        validate(need(!(is.null(input$DEGroup) || input$DEGroup == ""), ""))    
        vsd <- DESeq2::varianceStabilizingTransformation(dataDESeq())
        datapca <- DESeq2::plotPCA(vsd, intgroup=input$DEGroup, returnData = TRUE)
        percentVar <- round(100 * attr(datapca, "percentVar"))
        p <- ggplot2::ggplot(datapca, ggplot2::aes_string(x="PC1", y="PC2", color=input$DEGroup)) + 
                ggplot2::geom_point(aes(text=rownames(datapca)), size=1) + 
                ggplot2::xlab(paste0("PC1: ", percentVar[1], "% variance")) +
                ggplot2::ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
                ggplot2::ggtitle("Principal Component Analysis") +
                ggplot2::theme(
                    axis.text = element_text(size = 14),
                    axis.title = element_text(size = 14),
                    #legend.text = element_text(size = 12),
                    #legend.title = element_text(size = 16 ,face="bold"),
                    legend.title = element_blank(),
                    plot.title = element_text(size = 16, face="bold", hjust = 0.5),
                    aspect.ratio = 1
                )
        p <- plotly::ggplotly(p, type = "scatter")
        return(p)
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
        p <- degTab[, plotly::plot_ly(x = log2FoldChange, y = BHpvalue, text = rn, mode = "markers", color = group, type = "scatter")] 
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