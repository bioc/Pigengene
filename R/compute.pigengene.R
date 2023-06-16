compute.pigengene <- function(Data, Labels, modules, saveFile="pigengene.RData",
                              selectedModules="All", amplification=5,
                              doPlot=TRUE, verbose=0, dOrderByW=TRUE, naTolerance=0.05,
                              doWgcna=FALSE){
    ##
    ## modules: A vector of integers determining module assignments.
    ##^Named by column names of Data.
    message.if(me="Pigengenes...", verbose=verbose)
    pigengene <- list()
    ## Constants:
    varEpsilon <- 10^(-8)
    pigengene[["call"]] <- match.call()
    if ("All" %in% selectedModules) {
        selectedModules <- unique(modules)
    }
    if(is.null(names(Labels))){names(Labels) <- rownames(Data)}
    ## Files:
    if(!is.null(saveFile)){
        pBasename <- gsub(basename(saveFile), pattern="\\.RData$", replacement="_pvalue")
        pvalueCsvFile <- combinedPath(dir=dirname(saveFile), fn=paste(pBasename, ".csv", sep=""))
        membershipCsvFile <- combinedPath(dir=dirname(saveFile), fn=paste("membership.csv", sep=""))
        pigengene[["weightsCsvFile"]] <- membershipCsvFile
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
    pigengene[['Reptimes']] <- balanced$Reptimes
    origIds <- balanced$origSampleInds
    pigengene[["Data"]] <- Data
    rm(balanced, Data) ## Memory usage matters for big Data
    gc()
    ## Do all modules have at least one gene in columns of Data?
    hasNoGene <- !selectedModules %in% unique(modules[genes])
    if(any(hasNoGene)){
        stop(paste("There is no gene in columns of Data for the following modules.",
                   "Has selectedModules been set properly?!",
                   paste(selectedModules[hasNoGene], collapse=", ")))
    }
    ## Eigengenes:
    m1 <- paste("Computing eigengenes using", length(genes), "genes &", nrow(balancedData),
                "samples...")
    message.if(me=m1, verbose=verbose-1)
    ## Main computation:
    if(doWgcna){
        mERes <- WGCNA::moduleEigengenes(balancedData[ , genes, drop=FALSE],
                                         modules[genes], verbose=verbose-4,
                                         scale=TRUE)
        names(mERes$varExplained) <- colnames(mERes$eigengenes)
        rownames(mERes$eigengenes) <- rownames(balancedData)
        ## What if a module has only 1 gene? WGCNA return NaNs.
        sgms <- names(which(table(modules)==1)) ## Single gene modules
        sgms <- intersect(sgms, selectedModules)
        if(length(sgms)>0){
            sgmGenes <- genes[match(sgms, modules[genes])]
            mERes$eigengenes[, paste0("ME", sgms)] <- balancedData[, sgmGenes, drop=FALSE]
        }
        mERes$eigengenes <- as.matrix(mERes$eigengenes)
        pigengene[["eigenResults"]] <- mERes
        eigengenes <- as.matrix(mERes$eigengenes)
    } else {
        mERes <- list()
        for(m2 in unique(modules[genes])){
            genes2 <- names(which(modules==m2))
            ## Remove constant genes: 
            genes3 <- genes2[colSds(balancedData[ , genes2, drop=FALSE]) > varEpsilon]
            ## QC:
            if(length(genes3)==0)
                stop(paste("No variable gene in module", m2))
            message.if(paste("Running PCA for module", m2, ", which has",
                              length(genes3),"genes with some variation..."),
                       verbose=verbose-2)
            prcompRes <- prcomp(scale(balancedData[ , genes3, drop=FALSE]), scale.=TRUE)
            ##^ scale() is needed to avoid rare non-convergence by LAPACK,
            ##although scaling is already done at the top of prcomp()!
            ##https://stat.ethz.ch/pipermail/r-sig-ecology/2013-January/003493.html
            var2 <- (prcompRes$sdev[1])^2/sum(prcompRes$sdev^2)
            name2 <- paste0("ME", m2)
            mERes[["varExplained"]][name2] <- var2
            eigengene2 <- prcompRes$x[,"PC1"]
            eigengene2 <- eigengene2/sqrt(sum(eigengene2^2))
            message.if("Aligning module eigengene with average expression...",
                       verbose=verbose-2)
            averExpr <- rowMeans(balancedData[ , genes3, drop=FALSE], na.rm=TRUE)
            corAve  <- cor(averExpr, eigengene2, use="p")
            if(!is.finite(corAve))
                corAve <- 0
            if(corAve<0)
                eigengene2 <- -eigengene2
            mERes[["eigengenes"]] <- cbind(mERes$eigengenes, eigengene2)
            colnames(mERes[["eigengenes"]])[ncol(mERes[["eigengenes"]])] <- name2
            ## If a module has a single gene,
            ## the corresponding eigengene is the normalized expression of that genes.
        }        
        mERes$eigengenes <- as.matrix(mERes$eigengenes)
        pigengene[["eigenResults"]] <- mERes
        eigengenes <- as.matrix(mERes$eigengenes)
    }
    if(any(is.na(eigengenes)))
        warning("NA values in eigengenes!")
    message.if("Computing memberships...", verbose=verbose-1)
    n1 <- nrow(eigengenes)
    eigengenesBd <- as.matrix(eigengenes[rownames(balancedData), , drop=FALSE])
    membership <- suppressWarnings(stats::cor(balancedData, eigengenesBd))
    ## If a gene is almost constant, the above correlation may be NA although
    ## that gene is definitely in module 0,
    ## and NOT in any other module. Put it where it belongs to.
    membership[colSds(balancedData) < varEpsilon, ] <- as.numeric(colnames(membership)=="ME0")
    eigengenes <- eigengenes[origIds, , drop=FALSE]
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
        pigengene[["pvalues"]] <- pvalues
        pigengene[["log.pvalue"]] <- log.pvalue ## base 10
    }
    pigengene[["Labels"]] <- Labels
    pigengene[["eigengenes"]] <- eigengenesOrdered
    pigengene[["membership"]] <- membership
    pigengene[["orderedModules"]] <- modules[order(modules)]
    pigengene[["annotation"]] <- ann1
    pigengene[["saveFile"]] <- saveFile
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
        pigengene[["heavyToLow"]] <- rownames(membershipCsv)
    }
    if(!is.null(saveFile)){
        write.csv(file=membershipCsvFile, membershipCsv)
    }
    pigengene[["weights"]] <- membershipCsv

    ## The Pigengene output:
    class(pigengene) <- "pigengene"
    save.if(pigengene, file=saveFile, verbose=verbose-1)

    ## Plot:
    if(doPlot){
        if(is.null(saveFile)){
            stop("Although doPlot is TRUE, I cannot save the plots because saveFile is NULL !")
        }
        sf <- file.path(dirname(saveFile),"plots")
        plot.pigengene(x=pigengene, saveDir=sf,
                       selectedModules=selectedModules, verbose=verbose,
                       DiseaseColors="Auto")
    }
    return(pigengene)
}
