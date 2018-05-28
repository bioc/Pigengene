project.eigen <- function(
    Data, saveFile=NULL, pigengene, 
    naTolerance=0.05, verbose=0, ignoreModules=c())
{
    ## Projects expression on the given (MILE) eigenspace.
    ## Data: expression matrix with genes on columns.
    ## eigenFile: produced by compute.pigengene().
    ## naTolerance: In the range [0, 1]. If a gene (column of Data) has more than this
    ##^ ratio NAs, the gene will be ignored, otherwise the NA will be replaced by the
    ## mean expression of the gene.
    ## ignoreModules: these modules will be ignored, 
    ## e.g., if membership for module "ME-1" is not defined
    ##^set ignoreModules=c(-1).
    message.if("Projecting...", verbose=verbose)
    res <- list()
    Data <- as.matrix(Data)
    modules <- pigengene$orderedModules ## numeric, named by genes
    membership <- pigengene$membership ## (number of genes) * (number of modules)
    geoEigengenes <- pigengene$eigengenes
    p1 <- c() ## projected matrix
    notMatched <- c()
    tooNaGenes <- c()
    replacedNaNum <- 0
    for (m1 in unique(modules)){
        if(m1 %in% ignoreModules)
            next
        genes <- names(which(modules==m1))
        moduleName <- paste("ME", m1, sep="")
        matching <- match(genes, colnames(Data))
        notMatched <- c(notMatched, genes[is.na(matching)])
        genesMatched <- genes[!is.na(matching)]
        DataM <- Data[, genesMatched, drop=FALSE] ## Matched Data.
        ## Check for NA:
        checked <- check.nas(Data=DataM, naTolerance=naTolerance)  
        DataM <- checked$cleaned
        tooNaGenes <- c(tooNaGenes, checked$tooNaGenes)
        replacedNaNum <- sum(replacedNaNum, checked$replacedNaNum)
        ## Scaling:
        dataM <- scale(DataM)
        dataM[, colSds(DataM) < 10^(-8)] <- 0
        ## Inner product:
        if(! moduleName %in% colnames(membership)){
            stop(paste("membership for module", moduleName, "is not defined!"))
        }
        pm1 <- dataM %*% membership[colnames(dataM), moduleName, drop=FALSE]
        ## Normalizing:
        pm2 <- pm1 * (sqrt(sum(geoEigengenes[, moduleName]^2))/sqrt(sum(pm1^2)))
        p1 <- cbind(p1, pm2)
        colnames(p1)[ncol(p1)] <- moduleName
        if(any(is.na(matching)) & verbose >1)
            message(paste("In module", moduleName, sum(is.na(matching)), 
                        "genes out of", length(genes), "did not match."))
    }##End for (m1 in unique(modules)).
    ## All matched?
    notMatchedNum <- length(notMatched)
    if (notMatchedNum > 0 & verbose >0) 
        warning(paste(notMatchedNum, "genes could not match."))
    res[["tooNaGenes"]] <- tooNaGenes
    res[["replacedNaNum"]] <- replacedNaNum
    res[["projected"]] <- p1
    res[["notMatched"]] <- notMatched
    projected <- res
    if (!is.null(saveFile)) {
        save.if(projected, file=saveFile, verbose=verbose)
    }
    invisible(res)
}
