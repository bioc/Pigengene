draw.cor <- function(
    Data, savePath=".", doPlotDistCor=FALSE, doSaveCor=FALSE, 
    verbose=0)
{
    ##
    message.if("Correlations...", verbose=verbose)
    res <- list()
    dir.create(path=savePath, recursive=TRUE, showWarnings=FALSE)
    saveFile <- file.path(savePath, "plottedCor.RData")
    plotFile <- file.path(savePath, "cor.png")
    pCor <- cor(Data)
    inds <- which(var(Data)[, 1]==0)
    pCor[inds, ] <- pCor[, inds] <- 0
    corHeat <- pheatmap(abs(pCor))
    png(plotFile)
    pheatmap(pCor[corHeat$tree_row$order, corHeat$tree_col$order], 
             cluster_cols=FALSE, cluster_rows=FALSE)
    dev.off()
    res[["corHeat"]] <- corHeat
    if (doPlotDistCor) {
        dcorMatrix <- dcor.matrix(Data=Data)

        png(paste(savePath, "dcor.png", sep=""))
        dcorHeat <- pheatmap(dcorMatrix)
        dev.off()
        res[["dcorMatrix"]] <- dcorMatrix
        res[["dcorHeat"]] <- dcorHeat
        png(file.path(savePath, "corDcor.png"))
        pheatmap(pCor[dcorHeat$tree_row$order, dcorHeat$tree_col$order], 
                 cluster_cols=FALSE, cluster_rows=FALSE)
        dev.off()
    }
    res[["savePath"]] <- savePath
    res[["saveFile"]] <- saveFile
    res[["plotFile"]] <- plotFile
    if (doSaveCor) 
        res[["pCor"]] <- pCor
    plottedCor <- res
    save.if(plottedCor, file=saveFile, verbose=verbose)
    ##message.if(paste("plottedCor was saved in", saveFile), verbose=verbose)
    plottedCor[["pCor"]] <- pCor
    return(plottedCor)
}
