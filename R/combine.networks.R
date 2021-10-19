combine.networks <- function(nets, contributions, outPath, midfix="",powerVector=1:20,
                             verbose=1, RsquaredCut=0.75,  minModuleSize=5,  doRemoveTOM=TRUE,
                             datExpr, doReturNetworks=FALSE, doSave=FALSE, doIdentifyModule=TRUE){
                             

    message.if("Combining networks...", verbose=verbose)
    message.if(paste("RsquaredCut:", RsquaredCut), verbose=verbose-1)
    ## QC:
    if(length(contributions) != length(nets)){
        stop("nets and contributions must have the same length!")
    }
    result <- list()
    result[["call"]] <- match.call()
    result[["midfix"]] <- midfix
    nodes <- c()
    ## Compute the union of all nodes:
    for(ind in 1:length(nets)){
        nodes <- union(nodes, rownames(nets[[ind]]))
    }
    network <- matrix(0, nrow=length(nodes), ncol=length(nodes), dimnames=list(nodes, nodes))
    denominators <- matrix(0, nrow=length(nodes), ncol=length(nodes), dimnames=list(nodes, nodes))
    ## Combine the networks:
    for(ind in 1:length(nets)){
        netI <- nets[[ind]]
        network[rownames(netI), colnames(netI)] <- network[rownames(netI), colnames(netI)] +
            contributions[ind] * netI
        denominators[rownames(netI), colnames(netI)] <- denominators[rownames(netI), colnames(netI)] +
            contributions[ind]
    }
    network <- network/denominators

    ##Diagonal elements should be all 1
    diag(network) <- 1
                                  
    ## Recall identify modules function to return modules
    if (doIdentifyModule){
        identifiedMod <- identify.modules(network=network, outPath=outPath, midfix="",
                                       powerVector=powerVector, verbose=verbose-1,
                                       RsquaredCut=RsquaredCut, minModuleSize=minModuleSize,
                                       doRemoveTOM=doRemoveTOM, datExpr=datExpr, doSave=FALSE)
        result[["power"]] <- identifiedMod$power
        result[["fits"]] <- identifiedMod$fits
        result[["modules"]] <- identifiedMod$modules
        result[["net"]] <- identifiedMod$net
    }
    
    ## Clean up:
    if(doReturNetworks)
        result[["Network"]] <- network

    combinedNetwork <- result
    ##Save results
    if(doSave){
        file1 <- file.path(outPath, paste0("combinedNetwork", midfix, ".RData"))
        result[["combinedNetworkFile"]] <- file1
        save.if(x1=combinedNetwork, file=file1, verbose=verbose)
    }
    
    ## Output:
    return(combinedNetwork)
}
