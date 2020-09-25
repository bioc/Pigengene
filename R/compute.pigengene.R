compute.pigengene <- function(
    Data, Labels, modules, saveFile="pigengene.RData",
    selectedModules="All", amplification=5,
    doPlot=TRUE, verbose=0, dOrderByW=TRUE, naTolerance=0.05){
    ##
    ## modules: A vector of integers determining module assignments.
    ##^Named by column names of Data.
    message.if(me="Pigengenes...", verbose=verbose)
    result <- list()
    result[["call"]] <- match.call()
    if ("All" %in% selectedModules) {
        selectedModules <- unique(modules)
    }
    if(is.null(names(Labels))){names(Labels) <- rownames(Data)}
    ## Files:
    if(!is.null(saveFile)){
        pBasename <- gsub(basename(saveFile), pattern="\\.RData$", replacement="_pvalue")
        pvalueCsvFile <- combinedPath(dir=dirname(saveFile), fn=paste(pBasename, ".csv", sep=""))
        membershipCsvFile <- combinedPath(dir=dirname(saveFile), fn=paste("membership.csv", sep=""))
        result[["weightsCsvFile"]] <- membershipCsvFile
    }
    ## QC:
    c1 <- check.pigengene.input(Data=Data, Labels=Labels, na.rm=TRUE, naTolerance=naTolerance)
    Data <- c1$Data
    Labels <- c1$Labels
    Data <- Data[names(Labels), , drop=FALSE]
    ## Data:
    genes <- colnames(Data)[modules[colnames(Data)] %in% selectedModules]
    balanced <- balance(Data=Data, Labels=Labels,
                        amplification=amplification, verbose=verbose-1, naTolerance=naTolerance)
    balancedData <- balanced$balanced
    myDat <- balancedData[ , genes, drop=FALSE]
    result[['Reptimes']] <- balanced$Reptimes
    ## Do all modules have at least one gene in columns of Data?
    hasNoGene <- !selectedModules %in% unique(modules[colnames(myDat)])
    if(any(hasNoGene)){
        stop(paste("There is no gene in columns of Data for the following modules.",
                   "Has selectedModules been set properly?!",
                   paste(selectedModules[hasNoGene], collapse=", ")))
    }
    ## Eigengenes:
    m1 <- paste("Computing eigengenes using", ncol(myDat), "genes &", nrow(myDat),
                "samples...")
    message.if(me=m1, verbose=verbose-1)
    ## Main computation:
    mERes <- WGCNA::moduleEigengenes(myDat, modules[colnames(myDat)], verbose=verbose-4,
                                     scale=TRUE)
    names(mERes$varExplained) <- colnames(mERes$eigengenes)
    rownames(mERes$eigengenes) <- rownames(myDat)
    ## What if a module has only 1 gene? WGCNA return NaNs.
    sgms <- names(which(table(modules)==1)) ## Single gene modules
    sgms <- intersect(sgms, selectedModules)
    if(length(sgms)>0){
        mERes$eigengenes[, paste0("ME", sgms)] <- myDat[, match(sgms, modules[colnames(myDat)]), drop=FALSE]
    }
    mERes$eigengenes <- as.matrix(mERes$eigengenes)
    result[["eigenResults"]] <- mERes
    eigengenes <- as.matrix(mERes$eigengenes)
    if(any(is.na(eigengenes)))
        warning("NA values in eigengenes!")
    message.if("Computing memberships...", verbose=verbose-1)
    n1 <- nrow(eigengenes)
    membership <- stats::cor(balancedData,
                             as.matrix(eigengenes[rownames(balancedData), , drop=FALSE]))
    ## If a gene is almost constant, the above correlation may be NA altough
    ## that gene is definitely in module 0,
    ## and NOT in any other module.
    membership[colSds(balancedData) < 10^(-8), ] <- as.numeric(colnames(membership)=="ME0")
    eigengenes <- eigengenes[balanced$origSampleInds, , drop=FALSE]
    ann1 <- as.character(Labels[rownames(eigengenes)])
    ##^ pheatmap cannot work with e.g., TRUE
    ann1 <- as.data.frame(ann1)
    row.names(ann1) <- rownames(eigengenes)
    colnames(ann1) <- "Condition"
    matched <- match(paste("ME", selectedModules, sep=''), colnames(eigengenes))
    eigengenesOrdered <- eigengenes[, matched, drop=FALSE]

    ## Pvalues:
    if(length(unique(unlist(Labels)))>1){ ## more than one Label
        pvalues.function <- welch.pvalue ## Can be pvalues.manova or welch.pvalue
        message.if("Eigengene Pvalues ...", verbose=verbose-1)
        pvalues <- pvalues.function(Data=as.matrix(eigengenes), Labels[rownames(eigengenes)])
        log.pvalue <- as.data.frame(log10(as.numeric(pvalues$pvals[, "Bonferroni"])))
        row.names(log.pvalue) <- colnames(eigengenes)
        colnames(log.pvalue) <- "pvalue(log)"
        Size <- table(modules)
        names(Size) <- paste("ME", names(Size), sep="")
        pvalCsv <- cbind(Size[rownames(pvalues$pvals)], pvalues$pvals)
        colnames(pvalCsv)[1] <- "Size"
        if(!is.null(saveFile)){
            write.csv(pvalCsv, file=pvalueCsvFile)
        }
        result[["pvalues"]] <- pvalues
        result[["log.pvalue"]] <- log.pvalue ## base 10
    }
    result[["Data"]] <- Data
    result[["Labels"]] <- Labels
    result[["eigengenes"]] <- eigengenesOrdered
    result[["membership"]] <- membership
    result[["orderedModules"]] <- modules[order(modules)]
    result[["annotation"]] <- ann1
    result[["saveFile"]] <- saveFile
    membershipCsv <- cbind(membership, modules[rownames(membership)])
    colnames(membershipCsv)[ncol(membershipCsv)] <- "Module"
    Weight <- c() ## Will add this as a column to the CSV file.
    for(m1 in selectedModules){
        g1 <- names(which(modules==m1))
        g1 <- intersect(g1, rownames(membershipCsv))
        ##^ Maybe columns of Data are fewer than the length of modules.
        Weight[g1] <- membershipCsv[g1, paste0("ME",m1), drop=FALSE]
    }
    membershipCsv <- cbind(membershipCsv, "Weight"=Weight[rownames(membershipCsv)])
    if(dOrderByW){
        ordered <- c()
        toBeSortedModules <- intersect(selectedModules, membershipCsv[,"Module"])
        for(m1 in sort(unique(toBeSortedModules))){
            ## Consider the genes that are in module m1:
            membership1 <- membershipCsv[membershipCsv[,"Module"]==m1, , drop=FALSE]
            ## Compute the absolute weight of the above genes in module m1:
            absolute1 <- abs(membership1[ ,paste0("ME",m1)])
            ## Order based on absolute wight and add them at the bottom of the ordered matrix:
            ordered <- rbind(ordered, membership1[order(absolute1, decreasing=TRUE), ,drop=FALSE])
        }
        membershipCsv <- ordered
        result[["heavyToLow"]] <- rownames(membershipCsv)
    }
    if(!is.null(saveFile)){
        write.csv(file=membershipCsvFile, membershipCsv)
    }
    result[["weights"]] <- membershipCsv

    ## The Pigengene output:
    pigengene <- result
    class(pigengene) <- "pigengene"
    save.if(pigengene, file=saveFile, verbose=verbose-1)

    ## Plot:
    if(doPlot){
        if(is.null(saveFile)){
            stop("Although doPlot is TRUE, I cannot save the plots because saveFile is NULL !")
        }
        sf <- file.path(dirname(saveFile),"plots")
        ##dc <- c("red", "cyan", "green", "black", "pink", "brown", "yellow", "orange")
        ##dc <- dc[1:length(unique(Labels))]
        plot.pigengene(x=pigengene, saveDir=sf,
                       selectedModules=selectedModules, verbose=verbose,
                       DiseaseColors="Auto")
    }
    return(pigengene)
}
