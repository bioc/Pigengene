wgcna.one.step <- function(
    Data, power, saveDir=".", blockSize="All", saveTOMs=FALSE, doThreads=FALSE,
    verbose=0, seed=NULL)
{
    ##
    message.if(me="Identifying the modules (WGCNA)...", verbose=verbose)
    result <- list()
    result[["call"]] <- match.call()
    if (blockSize=="All") 
        blockSize <- ncol(Data)
    if(doThreads)
        WGCNA::allowWGCNAThreads()
    options(stringsAsFactors=FALSE)
    message.if(me=paste("power=", power), verbose=verbose)
    net <- blockwiseModules(datExpr=Data, power=power, TOMType="unsigned", 
                            minModuleSize=20, reassignThreshold=1e-06, 
                            mergeCutHeight=0.15, numericLabels=TRUE, 
                            pamRespectsDendro=FALSE, saveTOMs=saveTOMs, 
                            saveTOMFileBase=paste(saveDir, "DataTOM", sep=""), 
                            verbose=verbose-2, maxBlockSize=blockSize, 
                            randomSeed=seed)
    names(net$colors) <- colnames(Data)
    modules <- net$colors
    result[["net"]] <- net
    result[["genes"]] <- colnames(Data)
    moduleColors <- labels2colors(modules)
    result[["modules"]] <- modules
    result[["moduleColors"]] <- moduleColors
    result[["power"]] <- power
    if(verbose>0){
        message(paste(length(unique(modules)), "modules were identified:"))
        print(table(modules))
    }
    sizeFile <- NULL
    ## Saving:
    if(!is.null(saveDir)){
        plotFile <- file.path(saveDir, "plots", "dendro.png")
        dir.create(path=dirname(plotFile), recursive=TRUE, showWarnings=FALSE)
        netFile <- file.path(saveDir, "net.RData")
        save.if(net, file=netFile, verbose=verbose)
        result[["netFile"]] <- netFile
        ##Plot:
        ##sizeGrWindow(12, 9), an extra window left open.
        png(filename=plotFile)
        plotDendroAndColors(net$dendrograms[[1]],  addGuide=TRUE,
                            moduleColors[net$blockGenes[[1]]], 
                            "Module colors", dendroLabels=FALSE, hang=0.03, 
                            guideHang=0.05)
        dev.off()
        saveFile <- file.path(saveDir, "wgOneStep.RData")
        sizeFile <- file.path(saveDir,"plots","sizes.png")
        result[["sizeFile"]] <- sizeFile
        wgOneStep <- result
        save.if(wgOneStep, file=saveFile, verbose=verbose)
    }
    draw.module.sizes(plotFile=sizeFile, net=net, verbose=verbose)
    return(wgOneStep)
}
