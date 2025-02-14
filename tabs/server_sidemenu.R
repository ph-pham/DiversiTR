#-------------------------------------------------------------------------------------------------------------------------------------------#
#  Generate side menu section
#-------------------------------------------------------------------------------------------------------------------------------------------#  
# open if RepSeqDT() is loaded
observeEvent(is.RepSeqExperiment(RepSeqDT()), {
    # render menu show data
    output$showDataTab <- renderMenu({
    convertMenuItem(
        menuItem(tabName = "showDataTab",
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
                tabName = "computeDiversity")
        ), tabName = "showDataTab"
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
    output$multipleSampleTab <- renderMenu({
        convertMenuItem(
            menuItem(
                text = "comparative sample analysis",
                icon = icon("users", lib = "font-awesome"),
                tabName = "multipleSampleTab"
            ), 
            "multipleSampleTab"
        )
    })
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
    # show summary of the RepSeqExperiment object  
    output$summaryRDS <- output$summaryTXT <- renderUI({
        #flush.console()
        printHtml(RepSeqDT())
    })

    # library sizes
    output$histlibsizes <- renderPlot({
        cts <- RepSeq::assay(RepSeqDT()) 
        p1 <- histSums(cts[, sum(count), by=lib][,V1], xlab="Number of counts", ylab="Library size") + ggtitle("Library Sizes")
        p2 <- histSums(cts[, sum(count), by=VpJ][,V1], xlab="Number of counts", ylab="Number of clonotypes") + ggtitle("Clonotype Size")
        gridExtra::grid.arrange(p1, p2, ncol=2)
    })        
}) # end observeEvent RepSeqDT() is loaded
  
