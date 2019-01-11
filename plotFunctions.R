#------------------- plot functions --------------------#
#' plot a heatmap og V-J count 
#'
#' function plots a heatmap of V-J counts for sampleName 
#'
#' @param x a matrix of count data
#' @param sampleName sample to plot, sample name must be existed among x column names.
#' @return a heatmap of count
#' @export
#' @seealso \code{\link[pheatmap]{pheatmap}}
# @example
# plotCountVJ()
plotCountVJ <- function(x, sampleName=NULL, scale = c("counts", "percent", "cpm")) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.") 
    sNames <- rownames(sData(x))
    if (is.null(sampleName)) {
        index <- sNames[1]
        cat("Plot for the first sample in x:", index, ".\n")
    } else {
        index <- sampleName 
        if (length(sampleName) > 1) {
            index <- sampleName[1]
            cat("Only the first sample is used, sample:", index, ".\n")           
            }                  
        if (is.na(match(index, sNames))) {
            stop("Sample ", index, " not found in x.")                 
        }
    }
    tmp <- dcast(assay(x)[lib==index, ], J~V, value.var="count", fun.aggregate = sum)
    data2plot <- data.frame(tmp, row.names=1)
    
    if(scale == "percent"){
      data2plot <- prop.table(data2plot)
    }
    if(scale == "cpm"){
      data2plot <- prop.table(data2plot)*10^6
    }
   
    graph.title <- paste0("sample: ", index)
    # print(data2plot)
    if (requireNamespace("pheatmap", quietly = TRUE)) {
        p = pheatmap::pheatmap(data2plot, main=graph.title, cluster_rows = FALSE, cluster_cols = FALSE, cellheight = 20)#, cellwidth = 12, cellheight = 12)
        }
    return(p)
}

# plot usage
#
# plots a heatmap of percentages for each feature: VpJ/V/J/V-J of several repertoires.
#
# @param x a matrix of counts features in rows and samples in columns.
# @param graph.title a string . Default is "Segment usage". 
# @param cluster if TRUE, hierarchical clustering tree is show in the sample dimension using Euclidean distance and Ward's criterion for aggregation.
# @param ... others options will be passed through pheatmap
# @return a heatmap
# @seealso \code{\link[pheatmap]{pheatmap}}
#plotUsage <- function(x, graph.title="Segment usage", cluster=TRUE,...) {
#    if (missing(x)) stop("Dat set is missing")
#    if (!is.matrix(x)) x <- as.matrix(x)
#    if (cluster) {
#         pheatmap::pheatmap(x, main=graph.title, cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "ward.D", clustering_distance_cols = "euclidean",...)
#            } else {
#            pheatmap::pheatmap(x, main=graph.title, cluster_rows = FALSE, cluster_cols = FALSE, heigh = 3,...)
#        }       
#}

# plot a heatmap of diversity indices.
#
# plotDiversity plots a heatmap of V-J counts for sampleName .
#
# @param x a matrix of count data.
# @param graph.title a a title to graph.
# @param cluster if TRUE, hierarchical clustering tree is show in the sample dimension using Euclidean distance and Ward's criterion for aggregation.
# @param sample.annotation a data frame containing sample-related information, sample name must be existed among x column names.
# @return a heatmap of count
# @export
# @seealso \code{\link[pheatmap]{pheatmap}}
# @example
# plotCountVJ()
#plotDiversity <- function(x, graph.title="Shannon diversity", cluster=TRUE, sample.annotation=NA) {
#    if (missing(x)) stop("Diversity index data set is missing.")
#    if (!is.matrix(x)) x <- as.matrix(x)
#    if (!is.na(sample.annotation)) sample.annotation <- data.frame(sample.annotation) 
#    if (cluster) {
#    pheatmap::pheatmap(x, main=graph.title, 
#        cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "ward.D", clustering_distance_cols = "euclidean",
#        annotation_col=sample.annotation)
#    } else {
#        pheatmap::pheatmap(x, main=graph.title, cluster_rows = FALSE, cluster_cols = FALSE, annotation_col=sample.annotation)
#    } 
#}

#' plot spectratype
#'
#' function plots spectratype for V genes of a repertoire.
#'
#' @param x an object of class [\code{\linkS4class{RepSeqExperiment}}]
#' @param sampleName sample to plot, if NULL the first sample of x is plotted.
#' @return a barplot
#' @export
# @example
# plotSpectratype
plotSpectratype <- function(x, sampleName=NULL, scale = c("counts", "percent", "cpm")) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is epxected.") 
    CDR3length=pep=lib=CDR3aa <- NULL
    sNames <- rownames(sData(x))
    if (is.null(sampleName)) {
        index <- sNames[1]
        cat("Plot for the first sample in x:", index, ".\n")
    } else {
        index <- sampleName 
        if (length(sampleName) > 1) {
            index <- sampleName[1]
            cat("Only the first sample is used, sample:", index, ".\n")           
            }                  
        if (is.na(match(index, sNames))) {
            stop("Sample ", index, " not found in x.")                 
        }
    }
    # subset count
    scl <- match.arg(scale)
    if(scl == "counts"){
      data2plot <- assay(x)[lib == index, ][, CDR3length:=nchar(CDR3aa)][,.(.N), by = .(V, CDR3length)]#geom_bar
      p <- ggplot2::ggplot(data = data2plot, aes(x = CDR3length, y = N, fill = V))
    }
    if(scl == "percent"){
      data2plot <- assay(x)[lib == index, ][, CDR3length:=nchar(CDR3aa)][,.(.N), by = .(V, CDR3length)][,percent := prop.table(N)]#geom_bar
      p <- ggplot2::ggplot(data = data2plot, aes(x = CDR3length, y = percent, fill = V)) + scale_y_continuous(labels = scales::percent)
    }
    if(scl == "cpm"){
      data2plot <- assay(x)[lib == index, ][, CDR3length:=nchar(CDR3aa)][,.(.N), by = .(V, CDR3length)][,cpm := prop.table(N)*10^6]#geom_bar
      p <- ggplot(data = data2plot, aes(x = CDR3length, y = cpm, fill = V))
    }
    #ggplot
    p +
      geom_bar(stat = "identity", position="stack", colour = "black") +
      labs(title = paste(index,": V Distribution by aa length"), x = "CDR3 length (aa)", y = scl) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

#' plot Renyi's profiles
#'
#' function plots Renyi's entropy profile for all samples in a RepSeqExperiment object.
#'
#' @param x an object of class RepSeqExperiment
#' @param alpha a vector of parameters for estimating Renyi's entropy.  
#' @param level level of repertoire to estimated 
#' @param colorBy name of a factor in sample information data (run sData(x)) 
#' @return a graph
#' @export
# @example
# plotRenyiProfiles()
plotRenyiProfiles <- function(x, alpha=c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf), level=c("VpJ", "V", "J", "VJ", "CDR3aa"), colorBy=NULL) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    if (length(alpha) < 2) stop("At least 2 values of alpha is needed.")
    variable <- NULL
    sNames <- rownames(sData(x)) 
    levelChoice <- match.arg(level)
    # compute Renyi
    tmp <- renyiProfiles(x, scales=alpha, level=levelChoice)
    # plot 
    if (is.null(colorBy)) {
      #matplot
         #matplot(as.matrix(data2plot[, variable]), as.matrix(data2plot[, ..sNames]), type="l", 
         #xlab="alpha", ylab="Renyi's entropy", main=paste0("Level ", levelChoice,": Renyi profiles"))
      #ggplot
      colNames <- colnames(tmp)
      #setnames(data2plot, colnames(data2plot), c("alpha", sNames))
      data2plot <- melt(data = tmp, id.vars = "variable", measure.vars = sNames, variable.name = "lib")
      data2plot[, alpha := as.numeric(as.character(variable))]
      ggplot(data = data2plot, aes(x = alpha, y = as.numeric(value), color = lib, group = lib)) +
        geom_line() +
        labs(title = paste("Level ",levelChoice, ": Renyi's Entropy"), y = "Renyi's Entropy") +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
        } else {
            if (is.na(match(eval(colorBy), colnames(sData(x))))) stop(paste0(colorBy," not found in sData(x)."))
            #ggplot
            sdata <- sData(x)
            data2plot <- melt(data = tmp, id.vars = "variable", measure.vars = sNames, variable.name = "lib")
            data2plot[,paste(colorBy) := lapply(.SD, function(x){sdata[x, colorBy]}), .SDcols = "lib"]
            data2plot[, alpha := as.numeric(as.character(variable))]
            ggplot(data = data2plot, aes_string(x = "alpha", y = "value", colour = paste(colorBy))) +
              geom_line(aes(group = lib)) +
              labs(title = paste("Level ",levelChoice, ": Renyi's Entropy"), y = "Renyi's Entropy") +
              theme(plot.title = element_text(face = "bold", hjust = 0.5))

        }
}


#' plot V or J proportion
#'
#' function plots V or J proportion in all samples in a RepSeqExperiment object.
#'
#' @param x an object of class RepSeqExperiment
#' @param level level of repertoire to estimated 
#' @return a barplot
#' @export
# @example
# plotPropVJ

# plotPropVJ <- function(x , level = c("V", "J")){
#   if (missing(x)) stop("x is missing.")
#   if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
#   levelChoice <- match.arg(level)
#   p <- RepSeq::countFeatures(x, levelChoice)
#   pp <- melt(p, id = levelChoice, colnames(p)[2:4])
#   ggplot(data = pp, aes_string(x = levelChoice, y = "value", fill = "variable")) + geom_bar(stat = "identity", position="stack")
# }


plotPropVJ <- function(x , level = c("V", "J"), sample = NULL, scale = c("counts", "percent", "cpm")){
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  
  if (is.null(sample)) {
    sName <- rownames(sData(x))[1]
    cat("Plot for the first sample in x:", sName, ".\n")
  } else sName <- sample
  levelChoice <- match.arg(level)
  # p <- RepSeq::countFeatures(x, levelChoice)
  # p[,percent := .SD/sum(.SD), .SDcols = sName]
  
  data2plot <- data.table::copy(assay(x))[lib == sName, lapply(.SD, sum), by = levelChoice, .SDcols = "count"]
  if(scale == "percent"){
    data2plot[,percent := prop.table(count)]
    p <-   ggplot(data = data2plot, aes_string(x = levelChoice, y = "percent", fill = levelChoice)) +
      scale_y_continuous(labels = scales::percent)
  }
  if(scale == "cpm"){
    data2plot[,cpm := prop.table(count)*10^6]
    p <-   ggplot(data = data2plot, aes_string(x = levelChoice, y = "cpm", fill = levelChoice)) 
  }
  if(scale == "counts"){
    
    p <- ggplot(data = data2plot, aes_string(x = levelChoice, y = "count", fill = levelChoice)) 
  }
  p  +
    geom_bar(width = 0.7, stat = "identity", show.legend=F) +
    scale_x_discrete(limits = data2plot[order(-count)][[levelChoice]]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))+
    labs(title = paste(sample, ": ", levelChoice, " Distribution")) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

plotFreqVpJ <- function(x, sample){
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")

  if (is.null(sample)) {
    sName <- rownames(sData(x))[1]
    cat("Plot for the first sample in x:", sName, ".\n")
  } else sName <- sample
  data2plot <- data.table::copy(assay(x))
  #print(data2plot)
  f <- function(x){
    if(x==1) "1"
    else if(x<=10) "]1, 10]"
    else if(x<=100) "]10, 100]"
    else if(x<=1000) "]100, 1000]"
    else if(x<=10000) "]1000, 10000]"
    else "]10000, 100000]"
    }
   data2plot <- data2plot[lib == sName, lapply(.SD, sum), by = VpJ, .SDcols = "count"][,interval := unlist(lapply(count, f))]
   breaks <- unique(data2plot[,interval])
   plotBreaks <- breaks[order(nchar(breaks), breaks)]
   data2plot <- data2plot[,lapply(.SD, sum), by = interval, .SDcols = "count"][,percent := prop.table(count)]
  #data2plot[,interval] <-  factor(data2plot[,interval], levels = c("1", "]1, 10]", "]10, 100]", "]100, 1000]", ">1000"))
  #ggplot(data2plot, aes(x = count)) + geom_histogram(breaks = c(1, 10, 100, 1000, max(data2plot[,count]))) + scale_x_log10(breaks = c(1, 10, 100, 1000, max(data2plot[,count]))) + scale_y_log10() + labs(title = "TR sequence distribution", y = "# of occurences") + theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggplot(data = data2plot, aes(x = interval, y = percent, fill = interval)) +
    geom_bar(stat = "identity", show.legend=F) +
    scale_x_discrete(limits=plotBreaks) +
    scale_y_continuous(labels = scales::percent) +
    geom_text(aes(y=percent, label= paste(100*round(percent, 3), "%")), vjust=-1) +
    #coord_polar("y", start = 0) +
    labs(title = "TR sequence distribution", x = "count", y = "# of occurences") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

plotSpectratypebis <- function(x, sampleName=NULL, scale = c("counts", "percent", "cpm"), showCDR3 = F) {
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is epxected.") 
  CDR3length=pep=lib=CDR3aa <- NULL
  sNames <- rownames(sData(x))
  if (is.null(sampleName)) {
    index <- sNames[1]
    cat("Plot for the first sample in x:", index, ".\n")
  } else {
    index <- sampleName 
    if (length(sampleName) > 1) {
      index <- sampleName[1]
      cat("Only the first sample is used, sample:", index, ".\n")           
    }                  
    if (is.na(match(index, sNames))) {
      stop("Sample ", index, " not found in x.")                 
    }
  }
  scl <- match.arg(scale)
  # subset count
  if (scl == "counts"){
    tmp <- assay(x)[lib == index,][, CDR3length:=nchar(CDR3aa)]
    data2plot <- tmp[,.(.N), by = .(V, CDR3length)]
    p <- ggplot2::ggplot(data = data2plot, aes(x = CDR3length, y = N))
  }
  if (scl == "percent"){
    data2plot <- assay(x)[lib == index,][, CDR3length:=nchar(CDR3aa)][,.(.N), by = .(V, CDR3length)][,percent := prop.table(N), by = V]#geom_bar
    p <- ggplot2::ggplot(data = data2plot, aes(x = CDR3length, y = percent)) + scale_y_continuous(labels = scales::percent)
  }
  if (scl == "cpm"){
    data2plot <- assay(x)[lib == index,][, CDR3length:=nchar(CDR3aa)][,.(.N), by = .(V, CDR3length)][,cpm := prop.table(N)*10^6, by = V]#geom_bar
    p <- ggplot2::ggplot(data = data2plot, aes(x = CDR3length, y = cpm))
  }
  p <- p + ggplot2::geom_bar(stat = "identity") + scale_x_continuous(breaks = data2plot[, unique(CDR3length)]) +
        ggplot2::labs(title = paste(index,": V Distribution by aa length"), x = "CDR3 length (aa)", y = scl) +
        ggplot2::theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
        ggplot2::facet_wrap(~ factor(V, naturalsort::naturalsort(unique(V))), ncol = 4, scales = "free")
  print(p)
}


# plotSimilarityMatrix <- function(x, level=c("VpJ", "V", "J", "VJ", "CDR3aa"), colorBy = NULL) {
#   if (missing(x)) stop("x is missing.")
#   if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
#   variable <- NULL
#   levelChoice <- match.arg(level)
#   cols <- c("lib", levelChoice, "count")
#   tmp <- assay(x)[,..cols]
#   sdata <- sData(x)
#   if(!is.null(colorBy)){
#   if (is.na(match(eval(colorBy), colnames(sData(x))))) stop(paste0(colorBy," not found in sData(x)."))
#   sNames <- unique(sdata[, colorBy])
#   tmp[,paste(colorBy) := lapply(.SD, function(x){sdata[x, colorBy]}), .SDcols = "lib"]
#   data <- dcast(data = tmp, paste(levelChoice, "~", colorBy), value.var = "count")
#   matrixData <- as.matrix(data[,..sNames])
#   simmat <- proxy::simil(x = matrixData, by_rows = F, diag = T, upper = T)
#   mat <- as.matrix(simmat, diag = 0)
#   diag(mat) <-1
#   graph.title <- paste0(colorBy," similarity heatmap : ", levelChoice)
#   if (requireNamespace("pheatmap", quietly = TRUE)) {
#     p = pheatmap::pheatmap(mat, main=graph.title, cluster_rows = FALSE, cluster_cols = FALSE, cellheight = 15)#, cellwidth = 12, cellheight = 12)
#   }
#   }
#   else{
#     sNames <- rownames(sdata)
#     data <- dcast(data = tmp, paste(levelChoice, "~lib"), value.var = "count")
#     matrixData <- as.matrix(data[,..sNames])
#     simmat <- proxy::simil(x = matrixData, by_rows = F, diag = T, upper = T)
#     mat <- as.matrix(simmat, diag = 0)
#     diag(mat) <-1
#     graph.title <- paste(" similarity heatmap : ", levelChoice)
#     if (requireNamespace("pheatmap", quietly = TRUE)) {
#       p = pheatmap::pheatmap(mat, main=graph.title, cluster_rows = FALSE, cluster_cols = FALSE, cellheight = 15)#, cellwidth = 12, cellheight = 12)
#     }
#   }
#   return(p)
# }
plotDissimilarityMatrix <- function(x, level=c("VpJ", "V", "J", "VJ", "CDR3aa"), method = c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis"), binary="FALSE") {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    variable <- NULL
    levelChoice <- match.arg(level)
    methodChoice <- match.arg(method)
    cols <- c("lib", levelChoice, "count")
    tmp <- assay(x)[,..cols]
    sdata <- sData(x)
    sNames <- rownames(sdata)
    groups <- sdata[,unlist(lapply(sdata, is.factor)), drop = F]
    data <- dcast(data = tmp, paste(levelChoice, "~lib"), value.var = "count", fun.aggregate = sum)
    simmat <- data[, vegan::vegdist(t(.SD), method=methodChoice, diag=TRUE, upper=TRUE, binary=as.logical(binary)), .SDcols=sNames]
    graph.title <- paste0("dissimilarity heatmap : ", levelChoice)
    if (requireNamespace("pheatmap", method = methodChoice, quietly = TRUE)) {
      p <- pheatmap::pheatmap(simmat, main=graph.title, cluster_rows = TRUE, cluster_cols = TRUE, 
        treeheight_row = 0L, cellheight = 7, cellwidth = 7,
        annotation_col=groups[-1], annotation_row = groups[-1], show_rownames=FALSE, clustering_method = "ward.D")
      }
  return(p)
}

# plotFrequencySpectrum <- function(x, groupBy = FALSE, colorBy = NULL){
#   if (missing(x)) stop("x is missing.")
#   if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
#   if(!is.null(colorBy)){
#     if (is.na(match(eval(colorBy), colnames(sData(x))))) stop(paste0(colorBy," not found in sData(x)."))
#     sdata <- sData(x)
#     counts <- assay(x)
#     if(groupBy){
#       counts[,paste(colorBy) := lapply(.SD, function(x){sdata[x, colorBy]}), .SDcols = "lib"]
#       counts <- counts[,lapply(.SD, sum), by = c("VpJ", paste(colorBy)), .SDcols = "count"][,.N, by = c(paste(colorBy), "count")]
#       print(head(counts))
#       ggplot(data = counts, aes_string(x = "count", y = "N", colour = paste(colorBy))) + geom_point() +geom_line()+ scale_x_log10()
#     }else{
#       counts <- counts[,lapply(.SD, sum), by = c("VpJ", "lib"), .SDcols = "count"]
#       counts[,paste(colorBy) := lapply(.SD, function(x){sdata[x, colorBy]}), .SDcols = "lib"]
#       data2plot <- counts[,.N, by = c("lib", "count", paste(colorBy))]
#       return(ggplot(data = data2plot, aes_string(x = "count", y = "N", colour = paste(colorBy))) + geom_point() +geom_line(aes(group = lib))+ scale_x_log10())
# 
#     }
#   }
#   else{
#     counts <- assay(x)[,lapply(.SD, sum), by = .(VpJ, lib), .SDcols = "count"][,.N, by = .(lib, count)]
#     ggplot(data = counts, aes(x = count, y = N, colour = lib)) + geom_point() +geom_line()+ scale_x_log10()
#   }
# }


plotFrequencySpectrum <- function(x, groupBy = FALSE, colorBy = NULL){
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  counts <- assay(x)[,lapply(.SD, sum), by = .(VpJ, lib), .SDcols = "count"][,.N, by = .(lib, count)]
  if(!is.null(colorBy)){
    if (is.na(match(eval(colorBy), colnames(sData(x))))) stop(paste0(colorBy," not found in sData(x)."))
    sdata <- sData(x)
    counts[,group := lapply(.SD, function(x){sdata[x, colorBy]}), .SDcols = "lib"]
    if(groupBy){
      counts[,sd := lapply(.SD, sd), by = c("count", "group"), .SDcols = "N"]
      counts <- counts[,lapply(.SD, mean), by = c("group", "sd", "count"), .SDcols = "N"]
      print(counts)
      ggplot(data = counts, aes(x = count, y = N, colour = group)) + geom_errorbar(aes(ymin=N-sd, ymax=N+sd), width=.1, position = "dodge") + geom_point() +geom_line()+ scale_x_log10()
    }else{
      print(counts)
      ggplot(data = counts, aes(x = count, y = N, colour = group)) + geom_point() +geom_line(aes(group = lib))+ scale_x_log10()
    }
  }
  else{
    ggplot(data = counts, aes(x = count, y = N, colour = lib)) + geom_point() +geom_line()+ scale_x_log10()
  }
}

plotRarefaction<- function(x, colorBy = NULL, groupBy = FALSE){#, step){
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  sdata <- sData(x)
  print(colorBy)
  sNames <- unique(sdata[, colorBy])
  if(is.null(colorBy)){
    counts <- assay(x)
    dcounts <- dcast(data = counts, VpJ,~lib, value.var = "count", fun.aggregate = sum)
    cumulcounts <- dcounts[, lapply(.SD, cumsum), .SDcols = colnames(dcounts[,!"VpJ"])]
  }
  else{
    if (is.na(match(eval(colorBy), colnames(sData(x))))) stop(paste0(colorBy," not found in sData(x)."))
    if(groupBy){
      counts <- assay(x)[,paste(colorBy) := lapply(.SD, function(x){sdata[x, colorBy]}), .SDcols = "lib"]
      dcounts <- dcast(data = counts, paste("VpJ~",colorBy), value.var = "count", fun.aggregate = sum)
      print(dcounts[,!"VpJ"])
      #cumulcounts <- dcounts[, lapply(.SD, cumsum), .SDcols = colnames(dcounts[,!"VpJ"])]
    }
    else{
      counts <- assay(x)
      dcounts <- dcast(data = counts, VpJ~lib, value.var = "count", fun.aggregate = sum)
      #cumulcounts <- dcounts[, lapply(.SD, cumsum), .SDcols = colnames(dcounts[,!"VpJ"])]
    }
  }
  #print(nbcols(counts[, lapply(.SD, function(x){x==0}), .SDcols =  colnames(dcounts[,!"VpJ"])]))
  # VpJlist <- as.vector(counts[, VpJ])
  #print(dcounts)
  rare <- rarefaction(dcounts[,!"VpJ"], 10)#, step)
  nblines <- ncol(rare)
  raredt <- data.table(rare)
  #raredt[, size := c(1, step*1:(nrow(raredt)-1))]
  rarem <- melt(raredt, measure.vars = colnames(raredt[,!"size"]))
  if(!is.null(colorBy) & !groupBy){
    rarem[,paste(colorBy) := lapply(.SD, function(x){sdata[x, colorBy]}), .SDcols = "variable"]
    return(ggplot(data = rarem[value!=0], aes_string(x = "size", y = "value", colour = colorBy)) + geom_line())
  }
  ggplot(data = rarem[value!=0], aes(x = size, y = value, colour = variable)) + geom_line()
  # tcounts = t(as.matrix(dcounts[,!"VpJ"]))
  # #print(head(tcounts))
  # raremax <- min(rowSums(tcounts))
  # rarecurve(tcounts, step = 10000, sample = raremax,
  #            label = FALSE)
}


plotDistribVpJ <- function(x, colorBy = NULL){
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if(is.null(colorBy)) stop("need group for now")
  
  sdata <- sData(x)
  counts <- assay(x)
  counts[,group := lapply(.SD, function(x){sdata[x, colorBy]}), .SDcols = "lib"]
  print(counts)
  counts <- counts[, lapply(.SD, sum), by = .(group, VpJ), .SDcols = "count"]
  # counts[,order(-count), by = .(group, VpJ)]
  # #ggplot(data = counts, aes_string(x = ))
  # print(counts)

  #countsb <- frankv(counts, cols = "count", ties.method = "average")
  counts[,rank := lapply(.SD, frankv, ties.method = "min", order = -1L), by = group, .SDcols = "count"]
  counts <- unique(counts[,!"VpJ"])
  ggplot(data = counts, aes(x = rank, y = count, colour = group)) + geom_point() + scale_x_log10() + scale_y_log10()
}

plotVenn <- function(x, level = c("V", "J", "VJ", "VpJ", "CDR3aa"), colorBy = NULL){
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if(length(unique(sData(x)[,colorBy]))>5) stop("Too many groups or samples")
  levelChoice = match.arg(level)
  sdata <- sData(x)
  sNames <- unique(sdata[, colorBy])
  counts <- assay(x)
  counts[,group := lapply(.SD, function(x){sdata[x, colorBy]}), .SDcols = "lib"]
  dcounts <- dcast(data = counts, paste(levelChoice,"~group"), value.var = "count", fun.aggregate = sum)
  a <- limma::vennCounts(dcounts[,..sNames])
  limma::vennDiagram(a, circle.col = c("red", "blue", "green3", "purple", "yellow"))
  
}

plotmuScore <- function(x, level=c("V", "J", "VJ", "VpJ", "CDR3aa"), type=c("count", "usage")) {
  levelChoice = match.arg(level)
  typeChoice = match.arg(type)
  sdata <- sData(x)
  groups <- sdata[,unlist(lapply(sdata, is.factor)), drop = F]
  temp <- muScore(x, levelChoice, typeChoice)
  data2plot <- data.frame(temp, row.names=1)
  pheatmap::pheatmap(data2plot[, -ncol(data2plot)], cluster_rows=FALSE, cluster_cols = TRUE, annotation_col=groups[-1], show_rownames=TRUE, cellheight = 12,
                     color = colorRampPalette(c("lightgrey", "red"))(100))
}

plot2v2count <- function(x, level=c("V", "J", "VJ", "VpJ", "CDR3aa"), libs = NULL, scale = c("counts", "log")){
  levelChoice = match.arg(level)
  scaleChoice = match.arg(scale)
  cols <- c("lib", levelChoice, "count")
  counts <- assay(x)[lib ==libs[1] | lib == libs[2], ..cols]
  data2plot <- dcast(counts, paste(levelChoice, "~lib"), value.var="count", fun.aggregate = sum)
  p <- ggplot(data2plot, aes_string(libs[1], libs[2])) +
    geom_count(aes(alpha = ..n..),
               size = 2) + scale_alpha(range = c(0.25, 1)) +
    labs(title = paste("counts comparison per ", levelChoice)) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
  if(scale == "log")
    return(p + scale_x_log10() + scale_y_log10() + annotation_logticks())
  return(p)
    
}

brush2v2count <- function(x, level=c("V", "J", "VJ", "VpJ", "CDR3aa"), libs = NULL, plot = NULL){
  levelChoice = match.arg(level)
  counts <- assay(x)[lib ==libs[1] | lib == libs[2]]
  p <- dcast(counts, paste(levelChoice, "~lib"), value.var="count", fun.aggregate = sum)
  bp <- brushedPoints(p, plot)
  counts[eval(parse(text = levelChoice)) %in% bp[[levelChoice]]]
}