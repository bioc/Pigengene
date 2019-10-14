plot.pigengene <- function(x, saveDir=NULL, DiseaseColors=c("red", "cyan"),
                           fontsize=35, doShowColnames=TRUE, fontsizeCol=25,
                           doClusterCols=ncol(pigengene$eigengenes)>1, verbose=2,
                           doShowRownames="Auto",
                           pngfactor=max(2, ncol(pigengene$eigengenes)/16),
                           do0Mem=FALSE,
                           ...){
    result <- list()
    pigengene <- x
    ##QC:
    if(class(pigengene)!="pigengene")
        stop("pigengene must be of class 'pigengene' !")
    ## FILES:
    if(!is.null(saveDir)){
        dir.create(path=saveDir, recursive=TRUE, showWarnings=FALSE)
        message.if("Pigengene plots in:", verbose=verbose)
        message.if(saveDir, verbose=verbose)
        plotBaseName <- "pigengene"
        plotFile <- file.path(saveDir, paste(plotBaseName, ".png", sep=""))
        plotMemFile <- file.path(saveDir, paste("membership.png",sep=""))
        notRowsPlotFile <- gsub(plotFile, pattern="\\.png$", replacement="_notRows.png")
        pBasename <- paste(plotBaseName, "_pvalue", sep="")
        pvaluePlotFile <- file.path(saveDir, paste(pBasename, ".png", sep=""))
    }else{
        plotFile <- NULL
        plotMemFile <- NULL
        notRowsPlotFile <- NULL
        pvaluePlotFile <- NULL

    }
    conds <- unique(pigengene$annotation[, 1]) ##c(cond1, cond2)
    if(length(conds)!= length(DiseaseColors)){
        stop("The number of DiseaseColors must agree with the pigengene$annotation!")
    }
    names(DiseaseColors) <- conds
    log.pvalue <- pigengene$log.pvalue
    orderedModules <- pigengene$orderedModules
    membership <- pigengene$membership
    typeColor <- list(DiseaseColors)
    names(typeColor) <- colnames(pigengene$annotation)
    png2 <- function(aFile, pngf=pngfactor, ...){
        if(!is.null(aFile)){
            ## 2 times larger than normal png.
            png(aFile, width=480*pngf, height=480*pngf, ...)
        }else{
            dev.new()
        }
    }
    dof <- function(){
        if(!is.null(saveDir)){
            dev.off()
        }
    }
    if(doShowRownames=="Auto")
        doShowRownames <- nrow(pigengene$eigengenes) < 100
    ## Pvalues:
    pvaluePlot <- NA
    ## Don't add them to the heatmap unless
    if(!is.null(log.pvalue)){
        if(any(!is.finite(log.pvalue[, "pvalue(log)"])))
            stop("Not finite element in log.pvalue!")
        if(sd(rep(log.pvalue[, "pvalue(log)"], times=2))>0.001){
            pvaluePlot <- log.pvalue
            png2(aFile=pvaluePlotFile, pngf=1)
            plot(sort(as.numeric(log.pvalue[, 1])),
                 ylab="log10 p-value, Bonferroni",xlab="Eigengenes")
            dof()
            message.if(paste("log.pvalue was plotted in:", pvaluePlotFile),
                       verbose=verbose-2)
        }
    }
    ##
    png2(aFile=plotFile)
    heat <- pheatmap(pigengene$eigengenes, cluster_rows=TRUE,
                     cluster_cols=doClusterCols, fontsize=fontsize,
                     annotation_row=pigengene$annotation,
                     annotation_col=pvaluePlot, show_rownames=doShowRownames,
                     annotation_colors=typeColor)
    dof()
    m4 <- paste("eigengenes were plotted in:", plotFile)
    message.if(m4, verbose=verbose-2)
    ##
    genes <- rownames(membership)
    genes <- genes[order(pigengene$orderedModules[genes])]
    if(!do0Mem){
        genes <- intersect(genes, names(orderedModules[orderedModules!=0]))
    }
    png2(aFile=plotMemFile)
    pheatmap(abs(membership[genes, ]), annotation_colors=typeColor,
             cluster_rows=FALSE, cluster_cols=FALSE, fontsize=fontsize,
             show_rownames=FALSE)
    dof()
    message.if(paste("membership was plotted in:", plotMemFile), verbose=verbose-2)
    ##
    png2(aFile=notRowsPlotFile)
    heatNotRows <- pheatmap.type(Data=pigengene$eigengenes, cluster_cols=doClusterCols,
                                 fontsize=fontsize, annRow=pigengene$annotation,
                                 annotation_col=pvaluePlot,
                                 show_rownames=doShowRownames,
                                 annotation_colors=typeColor,
                                 show_colnames=doShowColnames,
                                 fontsize_col=fontsizeCol, ...)
    dof()
    message.if(paste("eigengenes (rows partially clustered) were plotted in:",
                     notRowsPlotFile), verbose=verbose-2)
    ##
    result[["heat"]] <- heat
    result[["heatNotRows"]] <- heatNotRows
    invisible(result)
}##End if(doPlot).
