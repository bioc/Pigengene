pheatmap.type <- function(Data, annRow, type=colnames(annRow)[1], ...){
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
        stop("annRow has rows that are not present in Data!")
    ## Put all samples in the same condition together.
    annRow <- annRow[order(annRow[, 1]), , drop=FALSE]
    Data <- Data[rownames(annRow), , drop=FALSE]
    conditions <- unique(as.character(annRow[, 1]))
    ##relData <- c() ## The relevant data.
    pheatmapS <- list()
    o1 <- c()
    for(cond in conditions){
        condSamples <- rownames(annRow)[which(annRow==cond)]  
        ##relData <- rbind(relData, Data[condSamples, , drop=FALSE])
        pa <- pheatmap(Data[condSamples, , drop=FALSE], cluster_cols=FALSE, silent=TRUE)
        pheatmapS[[as.character(cond)]] <- pa
        o1 <- c(o1, pa$tree_row$order+length(o1))
    }
    ann1 <- annRow[o1, , drop=FALSE]
    pAll <- pheatmap(Data[o1, , drop=FALSE], annotation_row=ann1, cluster_rows=FALSE, ...)
    res[["pheatmapS"]] <- pheatmapS
    res[["pheat"]] <- pAll
    res[["ordering"]] <- o1
    res[["annRowAll"]] <- ann1
    invisible(res)
}
