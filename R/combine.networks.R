combine.networks <- function(nets, contributions, outPath, midfix="", powerVector=1:20,
                             verbose=9, RsquaredCut=0.75, minModuleSize=5, doRemoveTOM=TRUE,
                             datExpr){

    ## QC:
    if(length(contributions) != length(nets)){
        stop("nets and contributions must have the same length!")
    }
    result <- list()
    result[["call"]] <- match.call()
    result[["midfix"]] <- midfix
    nodes <- c()
    ## compute the union of all nodes:
    for(ind in 1:length(nets)){
        nodes <- union(nodes, rownames(nets[[ind]]))
    }

    network <- matrix(0, nrow=length(nodes), ncol=length(nodes), dimnames=list(nodes, nodes))
    denominators <- matrix(0, nrow=length(nodes), ncol=length(nodes), dimnames=list(nodes, nodes))

    ## Combine the networks:
    for(ind in 1:length(nets)){
        netI <- nets[[ind]]
        network[rownames(netI), colnames(netI)] <- network[rownames(netI), colnames(netI)] +
            contributions[[ind]] * netI
        denominators[rownames(netI), colnames(netI)] <- denominators[rownames(netI), colnames(netI)] +
          contributions[[ind]]
    }
    result[["Network"]] <- network/denominators
    result[["denominators"]] <- denominators

    ##Diagonal elements should be all 1
    diag(result$Network) <- 1

    ##block size for both power estimation and blockwiseModule
    maxBlockSize <- nrow(result[["Network"]])
    ##power
    pstPower <- WGCNA::pickSoftThreshold.fromSimilarity(similarity=result$Network,
                                                        powerVector=powerVector,
                                                        blockSize=maxBlockSize,
                                                        verbose=verbose-1, RsquaredCut=RsquaredCut)
    result[["power"]] <- pstPower$powerEstimate
    result[["fits"]] <- pstPower$fitIndices

    ##Tom module for blockwiseModule function
    TOM <- as.dist(WGCNA::TOMsimilarity(adjMat=result$Network^pstPower$powerEstimate,
                               TOMType="unsigned", TOMDenom="min", verbose=TRUE))

    ##save TOM-module to load again for blockwiseModule
    tomFile <- file.path(outPath, paste0("TomModule", midfix, ".RData"))
    save(TOM, file=tomFile)
    file.remove("./blockwiseTOM-block.1.RData")
    command1 <- paste("ln -s ", tomFile, "./blockwiseTOM-block.1.RData")
    if(verbose>2)
        print(command1)
    system(command1)

    ## The following function won't use "network" but it can't work without it so,
    ## we should pass it as an argument. --Hanie.

## Changed name of object to "net" instead of "modules" to match what is done in wgcna.one.step() in Pigengene --Meg

    net <- WGCNA::blockwiseModules(datExpr=datExpr, saveTOMs=FALSE, power=pstPower$powerEstimate,   checkMissingData = FALSE,
                                  loadTOM=TRUE, maxBlockSize=maxBlockSize,
                                  mergeCutHeight=0.15, reassignThreshold = 1e-06,
                                  networkType="unsigned",
                                  replaceMissingAdjacencies=FALSE, TOMType="unsigned", TOMDenom="min",
                                  minModuleSize=minModuleSize, numericLabels=TRUE, verbose=TRUE)
    modules <- net$colors
    file.remove("./blockwiseTOM-block.1.RData")
    if(verbose>1)
        print(table(modules))
    names(modules) <- rownames(result$Network)
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

    if(doRemoveTOM)
      remove(file=tomFile)

    ##save results
    file1 <- file.path(outPath, paste0("combinedNetwork", midfix, ".RData"))
    result[["combinedNetworkFile"]] <- file1
    combinedNetwork <- result
    save(combinedNetwork, file=file1)
    if(verbose>0)
        print(paste("combinedNetwork was saved in", file1))
    return(combinedNetwork)
}