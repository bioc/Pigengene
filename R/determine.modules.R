determine.modules <- function(network, outPath, midfix="", powerVector=1:20,
                             verbose=1, RsquaredCut=0.75, minModuleSize=5, doRemoveTOM=FALSE,
                             datExpr, doSave=FALSE){
    ## It identifies modules of the  network by using of WGCNA(). 
    ## Network: is a matrix which is builded in combine.network().
    ## Output: modules, the identified modules of the network.

    message.if("Identifying modules...", verbose=verbose)
    result <- list()
    result[["call"]] <- match.call()
    result[["midfix"]] <- midfix
    nodes <- c()

    
    ##block size for both power estimation and blockwiseModule
    maxBlockSize <- nrow(network)
    ##power
    if(verbose<2)
        sink(nullfile()) ## To suppress the table printed by WGCNA.
    pstPower <- WGCNA::pickSoftThreshold.fromSimilarity(similarity=network,
                                                        powerVector=powerVector,
                                                        blockSize=maxBlockSize,
                                                        verbose=verbose-1, RsquaredCut=RsquaredCut)
    sink()
    message.if(paste("power:", pstPower$powerEstimate), verbose=verbose-1)
    result[["power"]] <- pstPower$powerEstimate
    result[["fits"]] <- pstPower$fitIndices
    if(is.na(pstPower$powerEstimate))
        stop("Consider a lower value for RsquaredCut, power is NA!")
    
    ##Tom module for blockwiseModule function
    TOM <- as.dist(WGCNA::TOMsimilarity(adjMat=network^pstPower$powerEstimate,
                                        TOMType="unsigned", TOMDenom="min", verbose=verbose-2))

    ##save TOM-module to load again for blockwiseModule
    tomFile <- file.path(outPath, paste0("TomModule", midfix, ".RData"))
    save(TOM, file=tomFile)
    tempFile <- "./blockwiseTOM-block.1.RData"
    if(file.exists(tempFile)){
        file.remove(tempFile)
    }
    command1 <- paste("ln -s ", tomFile, tempFile)
    if(verbose>2)
        print(command1)
    system(command1)

    ## The following function won't use "network" but it can't work without it so,
    ## we should pass it as an argument. --Hanie.
    ## Changed name of object to "net" instead of "modules" to match what is done
    ##in wgcna.one.step() in Pigengene --Meg
    net <- WGCNA::blockwiseModules(datExpr=datExpr, saveTOMs=FALSE, power=pstPower$powerEstimate,
                                   checkMissingData = FALSE,
                                   loadTOM=TRUE, maxBlockSize=maxBlockSize,
                                   mergeCutHeight=0.15, reassignThreshold = 1e-06,
                                   networkType="unsigned",
                                   replaceMissingAdjacencies=FALSE, TOMType="unsigned", TOMDenom="min",
                                   minModuleSize=minModuleSize, numericLabels=TRUE, verbose=verbose-2)
    modules <- net$colors
    if(file.exists(tempFile)){
        file.remove(tempFile)
    }
    if(verbose>1)
        print(table(modules))
    names(modules) <- rownames(network)
    result[["modules"]] <- modules
    result[["net"]] <- net

    ## Plot:
    png(filename=file.path(outPath, paste0("Module_sizes", midfix, ".png")),
        height=2*480, width=2*480, res=2.5*72)
    plot(log10(table(modules)[-1]), col="blue", ylab="Number of genes", xlab="Module index",
         cex.lab=1.4, yaxt="n")
    axis(2, at=c(0, 1, 2, 3), labels=c(1, 10, 100, 1000))
    dev.off()
    if(verbose>1)
        print(paste("There are", table(modules)[1], "outliers (ME0)"))

    ## Clean up:
    if(doRemoveTOM)
        unlink(tomFile)
    return(result)
}
