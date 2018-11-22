one.step.pigengene <- function(
    Data, saveDir="Pigengene", 
    Labels, testD=NULL, testLabels=NULL, doBalance=TRUE, RsquaredCut=0.8,
    costRatio=1, toCompact=FALSE, bnNum=0, bnArgs=NULL, useMod0=FALSE, 
    mit="All", ## unique(Labels)[1], 
    verbose=0, doHeat=TRUE, seed=NULL)
{
    ## costRatio: Implemented only for 2 classes.
    ##^Determines how severe it is to misclassify a sample accross types.
    ##^E.g., if costRatio=2, misclassification of a sample of the 1st type is
    ##considered twice worse than misclassification of a sample of the 2nd type.
    results <- list()
    results[["call"]] <- match.call()
    m1 <- paste("Pigengene started analizing", nrow(Data), 
                "samples using", ncol(Data), "genes...")
    message.if(me=m1, verbose=verbose)
    if(verbose>1){
        print(table(Labels))
    }
    ## QC:
    c1 <- check.pigengene.input(Data=Data, Labels=Labels, na.rm=TRUE)
    Data <- c1$Data
    Labels <- c1$Labels
    if(!is.null(testD)){
        ct <- check.pigengene.input(Data=testD, Labels=testLabels, na.rm=TRUE)
        testD <- ct$Data
        testLabels <- ct$Labels
    }
    ## saveDir:
    if(length(grep(saveDir, pattern=" ")>0))
        stop("saveDir cannot have space!")
    dir.create(path=saveDir, recursive=TRUE, showWarnings=FALSE)

    ## Data for WGCNA:
    wData <- switch(mit, 
                    "All"=Data, 
                    Data[which(Labels %in% mit), ])
    ##stop('mit must be equal to "cond1", "cond2", or "Both"!'))
    ## WGCNA:
    if(doBalance)
        wData <- balance(Data=wData, Labels=Labels, verbose=verbose-1)$balanced
    calculateBetaRes <- calculate.beta(saveFile=NULL, RsquaredCut=RsquaredCut,
                                       Data=wData, verbose=verbose-1)
    results[["betaRes"]] <- calculateBetaRes
    if(is.na(calculateBetaRes[["power"]]))
        stop("Consider a lower value for RsquaredCut, power is NA!")
    wgRes <- wgcna.one.step(Data=wData, seed=seed, 
                            power=calculateBetaRes[["power"]], 
                            saveDir=saveDir, verbose=verbose-1)
    rm(wData)
    results[["moduleRes"]] <- wgRes
    ## Eigengenes:
    pigengene <- compute.pigengene(Data=Data, Labels=Labels, 
                                   modules=wgRes$net$colors, 
                                   saveFile=combinedPath(saveDir, 'pigengene.RData'), 
                                   doPlot='TRUE', verbose=verbose)
    results[["pigengene"]] <- pigengene
    ## Multiple conditions?
    if(length(unique(Labels))<2){
        message.if("A single condition in Labels. BN and trees are skipped.",
                   verbose=verbose)
        return(results)
    }
    selectedFeatures <- NULL
    ## BN:
    if(bnNum !=0){
        ## Arguments:
        bnArgs$bnPath <- combinedPath(saveDir, 'bn')
        bnArgs$doShuffle <- if(!is.null(bnArgs$doShuffle)) bnArgs$doShuffle else TRUE
        bnArgs$tasks <- if(!is.null(bnArgs$tasks)) bnArgs$tasks else "All"
        bnArgs <- c(bnArgs, list(pigengene=pigengene, bnNum=bnNum, verbose=verbose-1, seed=seed))
        ## Call
        learnt <- do.call(learn.bn, bnArgs)
        results[["leanrtBn"]] <- learnt  
        BN <- learnt$consensus1$BN  
        results[["BN"]] <- BN
        ##^ The fist threshould is used for selecting.
        selectedFeatures <- c()
        if(learnt$use.Disease)
            selectedFeatures <- c(selectedFeatures, children("Disease", x=BN))
        if(learnt$use.Effect)
            selectedFeatures <- c(selectedFeatures, parents("Effect", x=BN))
        selectedFeatures <- setdiff(selectedFeatures, c("Disease", "Effect"))
        if(length(selectedFeatures)==0){
            warning("The condition variable has no child. BN results will be ignored.")
            selectedFeatures <- NULL
        }
    }

    ## Trees:
    c5Path <- combinedPath(saveDir, 'C5Trees')
    dir.create(path=c5Path, recursive=TRUE, showWarnings=FALSE)
    c5treeRes <- make.decision.tree(pigengene=pigengene, Data=Data, 
                                    testD=testD, testL=testLabels, 
                                    selectedFeatures=selectedFeatures, saveDir=c5Path, 
                                    minPerLeaf=NULL, useMod0=useMod0, doHeat=doHeat, 
                                    costRatio=costRatio, verbose=verbose, 
                                    toCompact=toCompact)
    ##
    results[["selectedFeatures"]] <- selectedFeatures
    results[["c5treeRes"]] <- c5treeRes
    return(results)
}
