pheatmap.type <- function(Data, annRow, type=colnames(annRow)[1], doTranspose=FALSE, conditions="Auto", ...){
    ## annRow: A data frame with row names the same as row names of Data.
    ## type: The column name of annRow representing two or more conditions.
    ## This function first performs hierarchical clustering on samples
    ## (rows of Data) within each condition.
    ##^Then plots a heatmap without further clustering of rows. 
    ## ... are passed to pheatmap function.
    res <- list()
    annRow <- annRow[, type, drop=FALSE]
    ## QC:
    if(is.null(rownames(annRow)))
        stop("annrow must have row names!")
    if(any(! rownames(annRow) %in% rownames(Data)))
        stop("annRow has rows that are not present in rows of Data!")
    ## Put all samples in the same condition together.
    annRow <- annRow[order(annRow[, 1]), , drop=FALSE]
    samplesOriginalOrder <- rownames(Data)
    ##Data <- Data[rownames(annRow), , drop=FALSE]
    if(conditions[1]=="Auto")
        conditions <- unique(as.character(annRow[, 1]))
    if(any(!conditions %in% unique(as.character(annRow[, 1])))){
        warning("Some of the conditions are not in annRow.")
        conditions <- intersect(conditions, unique(as.character(annRow[, 1])))
    }
    pheatmapS <- list()
    dataPlot <- c()
    ann1 <- c()
    for(cond in conditions){
        condSamples <- rownames(annRow)[which(annRow==cond)]  
        pa <- pheatmap(Data[condSamples, , drop=FALSE], cluster_cols=FALSE, silent=TRUE)
        pheatmapS[[as.character(cond)]] <- pa
        o2 <- pa$tree_row$order
        dataPlot <- rbind(dataPlot, Data[condSamples[o2], , drop=FALSE])
        ann1 <- rbind(ann1, annRow[condSamples[o2], , drop=FALSE])
    }
    if(!doTranspose){
        pAll <- pheatmap(dataPlot, annotation_row=ann1, cluster_rows=FALSE, ...)
    } else { ## Transpose
        pAll <- pheatmap(t(dataPlot), annotation_col=ann1, cluster_cols=FALSE, ...)
    }
    res[["pheatmapS"]] <- pheatmapS
    res[["pheat"]] <- pAll
    res[["ordering"]] <- match(rownames(dataPlot), samplesOriginalOrder)
    res[["annRowAll"]] <- ann1
    invisible(res)
}
