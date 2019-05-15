draw.cor.cond <- function(Data, Labels, savePath=".", verbose=0, naTolerance=0.05){
    message.if("Correlations in conditions...", verbose=verbose)
    res <- list()
    c1 <- check.pigengene.input(Data=Data, Labels=Labels,na.rm=TRUE, naTolerance=naTolerance)
    Data <- c1$Data
    Labels <- c1$Labels

    for(l1 in unique(Labels)){
        savePath1 <- combinedPath(savePath, paste("cor_",l1,"", sep=""))
        cond1Plot <- draw.cor(Data=Data[names(which(Labels==l1)),],
                              savePath=savePath1, verbose=verbose-1)
        if(FALSE){ ## Not needed for now.
            png(gsub(cond1Plot$plotFile, pattern="\\.png",
                     replacement="-other.png"))
            pheatmap(cond2Plot$pCor[cond1Plot$corHeat$tree_row$order,
                                    cond1Plot$corHeat$tree_col$order],
                     cluster_cols=FALSE, cluster_rows=FALSE)
            dev.off()
        }
        res[[paste(l1,"Plot",sep="")]] <- cond1Plot
    }
    return(res)
}
