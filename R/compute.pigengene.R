compute.pigengene <- function(
    Data, Labels, modules, saveFile="pigengene.RData", 
    selectedModules="All", amplification=5, 
    doPlot=TRUE, verbose=0)
{
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
    pBasename <- gsub(basename(saveFile), pattern="\\.RData$", replacement="_pvalue")
    pvalueCsvFile <- combinedPath(dir=dirname(saveFile), fn=paste(pBasename, ".csv", sep=""))
    membershipCsvFile <- combinedPath(dir=dirname(saveFile), fn=paste("membership.csv", sep=""))
    ## QC:
    c1 <- check.pigengene.input(Data=Data, Labels=Labels, na.rm=TRUE)
    Data <- c1$Data
    Labels <- c1$Labels
    Data <- Data[names(Labels), ]
    ## Data:
    genes <- names(modules)[modules[colnames(Data)] %in% selectedModules]
    origData <- Data[, genes, drop=FALSE]
    balanced <- balance(Data=origData, Labels=Labels, 
                        amplification=amplification, verbose=verbose-1)
    myDat <- balanced$balanced
    ## Eigengenes:
    m1 <- paste("Computing eigengenes using", ncol(myDat), "genes &", nrow(myDat), 
                "samples...")
    message.if(me=m1, verbose=verbose-1)
    ## Main computation:
    browser()
    eigenResults <- moduleEigengenes(myDat, modules[colnames(myDat)], verbose=verbose-4)
    names(eigenResults$varExplained) <- colnames(eigenResults$eigengenes)
    rownames(eigenResults$eigengenes) <- rownames(myDat)
    result[["eigenResults"]] <- eigenResults
    eigengenes <- eigenResults$eigengenes
    message.if("Computing memberships...", verbose=verbose-1)
    n1 <- nrow(eigengenes)
    membership <- cor(myDat, as.matrix(eigengenes[rownames(myDat), , drop=FALSE]))
    eigengenes <- eigengenes[balanced$origSampleInds, , drop=FALSE]
    ann1 <- as.character(Labels[rownames(eigengenes)])
    ##^ pheatmap cannot work with e.g., TRUE
    ann1 <- as.data.frame(ann1)      
    row.names(ann1) <- rownames(eigengenes)
    colnames(ann1) <- "Condition"
    matched <- match(paste("ME", selectedModules, sep=''), colnames(eigengenes))
    eigengenesOrdered <- eigengenes[, matched, drop=FALSE]
    ## Pvalues:
    if(length(unique(Labels))>1){
        message.if("Eigengene Pvalues ...", verbose=verbose-1)
        pvalues <- pvalues.manova(Data=as.matrix(eigengenes), Labels[rownames(eigengenes)])   
        log.pvalue <- as.data.frame(log10(as.numeric(pvalues$pvals[, "Bonferroni"])))
        row.names(log.pvalue) <- colnames(eigengenes)
        colnames(log.pvalue) <- "pvalue(log)"
        Size <- table(modules)
        names(Size) <- paste("ME", names(Size), sep="")
        pvalCsv <- cbind(Size, pvalues$pvals) 
        write.csv(pvalues$pvals, file=pvalueCsvFile)
        result[["pvalues"]] <- pvalues
        result[["log.pvalue"]] <- log.pvalue ## base 10
    }
    result[["Data"]] <- Data
    result[["Labels"]] <- Labels
    result[["eigengenes"]] <- eigengenesOrdered
    result[["membership"]] <- membership
    result[["orderedModules"]] <- modules[order(modules)]
    result[["annotation"]] <- ann1
    result[["weightsCsvFile"]] <- membershipCsvFile
    result[["saveFile"]] <- saveFile
    pigengene <- result
    class(pigengene) <- "pigengene"
    save.if(pigengene, file=saveFile, verbose=verbose-1)
    membershipCsv <- cbind(membership, modules[rownames(membership)])
    colnames(membershipCsv)[ncol(membershipCsv)] <- "Module" 
    write.csv(file=membershipCsvFile, membershipCsv)
    if(doPlot){
        sf <- file.path(dirname(saveFile),"plots")
        ##sf <- gsub(saveFile, pattern="\\.RData$", replacement="")
        dc <- c("red", "cyan", "green", "black", "pink", "brown", "yellow", "orange")
        dc <- dc[1:length(unique(Labels))]
        plot.pigengene(x=pigengene, saveDir=sf, 
                       selectedModules=selectedModules, verbose=verbose,
                       DiseaseColors=dc)
    } 
    return(pigengene)
}
