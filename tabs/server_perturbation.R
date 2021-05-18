#-------------------------------------------------------------------------------------------------------------------------------------------#
#  perturbation section
#-------------------------------------------------------------------------------------------------------------------------------------------#  
    # render select group UI
    output$PertGroupUI <- renderUI({
        selectGroup("PertGroupSelected", RepSeqDT())
    })
    # render selection control group
    output$CtrlGroupUI <- renderUI({
        validate(need(!(is.null(input$PertGroupSelected) || input$PertGroupSelected ==""), ""))
        selectizeInput("CtrlGroup",
            "Select control group",
            choices = levels(sData(RepSeqDT())[, input$PertGroupSelected]),
            options = list(onInitialize = I('function() { this.setValue(""); }')),
            multiple = F
        )
    })
    # render selection distance
    output$PertDistUI <- renderUI({
        validate(need(!(is.null(input$PertGroupSelected) || input$PertGroupSelected ==""), ""))
        validate(need(!(is.null(input$CtrlGroup) || input$CtrlGroup == ""), ""))
        selectizeInput("pertDist",
            "Select distance",
            choices = list("manhattan", "euclidean"),
            options = list(onInitialize = I('function() { this.setValue(""); }')),
            multiple = F
        )
    })
    # create perturbation data 
    dataPert <- reactive({
        validate(need(!(is.null(input$PertGroupSelected) || input$PertGroupSelected ==""), ""))
        validate(need(!(is.null(input$CtrlGroup) || input$CtrlGroup == ""), ""))
        validate(need(!(is.null(input$pertDist) || input$pertDist == ""), ""))
        sampleinfo <- sData(RepSeqDT())
        ctrnames <- rownames(sampleinfo)[which(sampleinfo[, input$PertGroupSelected] %in% input$CtrlGroup)]
        pertscore <- perturbation(x = RepSeqDT(), ctrl.names = ctrnames, distance = input$pertDist)
        return(pertscore)
    })
    # pca
    output$pertPCA <- renderPlot({
        validate(need(!(is.null(input$PertGroupSelected) || input$PertGroupSelected ==""), ""))
        validate(need(!(is.null(input$CtrlGroup) || input$CtrlGroup == ""), ""))
        validate(need(!(is.null(input$pertDist) || input$pertDist == ""), ""))
        perturb <- dataPert()
        pca <- prcomp(t(perturb), scale = TRUE)
        varexp <- (pca$sdev^2 / sum(pca$sdev^2))*100
        scores <- as.data.frame(pca$x)
        scores <- merge(scores, sData(RepSeqDT()), by = 0)
        p <- ggplot(scores, aes_string(x="PC1", y="PC2")) +
                geom_point(shape = 21, size = 2.5, stroke=0.5, aes_string(fill = input$PertGroupSelected, col = input$PertGroupSelected)) +
                theme_bw() +
                xlab(paste("PC1 (",round(varexp[1], digit=2),"%)", sep = ""))+
                ylab(paste("PC2 (",round(varexp[2], digit=2),"%)", sep = ""))+
                ggtitle("Principal Component Analysis")+
                theme(
                    axis.text = element_text(size = 16),
                    axis.title = element_text(size = 16),
                    legend.text = element_text(size = 16),
                    #legend.title = element_text(size = 16 ,face="bold"),
                    legend.title = element_blank(),
                    plot.title = element_text(size = 18, face="bold", hjust = 0.5),
                    aspect.ratio = 1,
                    legend.position = "bottom"
                )
        #p <- plotly::ggplotly(p, type = "scatter") 
        return(p)
    })
    # heatmap
    output$pertHeatmap <- renderPlot({
        validate(need(!(is.null(input$PertGroupSelected) || input$PertGroupSelected ==""), ""))
        perturb <- dataPert()
        sampleinfo <- sData(RepSeqDT())[,unlist(lapply(sData(RepSeqDT()), is.factor))]
        pheatmap::pheatmap(perturb, main = paste0("perturbation against: ", input$PertGroupSelected, ":", input$CtrlGroup), 
                            scale = "none", cluster_rows = TRUE, cluster_cols = TRUE, annotation_col = sampleinfo, annotation_legend = FALSE,  
                            show_colnames = TRUE, show_rownames = TRUE, clustering_method = "ward.D", fontsize = 8)
    })
    # data output
    output$PertTab <- renderDataTable({
        return(datatable(dataPert(), options = list(scrollX=TRUE, dom = 'Bfrtip', pageLength = 10)) %>% formatRound(c(1: ncol(dataPert())), 2))
    })
    