one.step.pigengene <- function(
    Data, saveDir="Pigengene",
    Labels, testD=NULL, testLabels=NULL, doBalance=TRUE, RsquaredCut=0.8,
    costRatio=1, toCompact=FALSE, bnNum=0, bnArgs=NULL, useMod0=FALSE,
    mit="All", ## unique(Labels)[1],
    verbose=0, doHeat=TRUE, seed=NULL, dOrderByW=TRUE)
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
    dataNum <- 1
    if(class(Data) =="list"){
        dataNum <- length(Data)
    }
    ## saveDir:
    if(length(grep(saveDir, pattern=" ")>0))
        stop("saveDir cannot have space!")
    dir.create(path=saveDir, recursive=TRUE, showWarnings=FALSE)

    if(!is.null(testD)){
        ct <- check.pigengene.input(Data=testD, Labels=testLabels, na.rm=TRUE)
        testD <- ct$Data
        testLabels <- ct$Labels
    }

    nets <- list()
    cont <- list()
    checkedData <- list()    
    checkedLabels <- list()    
    for(ind in 1:dataNum){
        if(dataNum==1){
	    DataI  <- Data
            LabelsI <- Labels
        } else {
            DataI <- Data[[ind]]
            LabelsI <- Labels[[ind]]
            ## contribution from the given data set DataI
	    cont[[ind]] <- nrow(DataI)
        }
        ## QC: 
	c1 <- check.pigengene.input(Data=DataI, Labels=LabelsI, na.rm=TRUE)
        DataI <- c1$Data
        LabelsI <- c1$Labels
	checkedData[[ind]] <- DataI
	checkedLabels[[ind]] <- LabelsI
        ## Data for WGCNA:
        wData <- switch(mit,
                        "All"=DataI,
                        DataI[which(LabelsI %in% mit), ])
        ##stop('mit must be equal to "cond1", "cond2", or "Both"!'))
        ## WGCNA:
        if(doBalance)
            wData <- balance(Data=wData, Labels=LabelsI, verbose=verbose-1)$balanced
        if(dataNum==1){
            calculateBetaRes <- calculate.beta(saveFile=NULL, RsquaredCut=RsquaredCut,
                                               Data=wData, verbose=verbose-1)
            results[["betaRes"]] <- calculateBetaRes
            if(is.na(calculateBetaRes[["power"]]))
                stop("Consider a lower value for RsquaredCut, power is NA!")
            betaI <- calculateBetaRes[["power"]]
            wgRes <- wgcna.one.step(Data=wData, seed=seed,
                                    power=betaI,
                                    saveDir=saveDir, verbose=verbose-1)
        } else {
            nets[[ind]] <- stats::cor(wData)
        }
        rm(wData)
    }



    ## Now use combine.networks()
    if(dataNum!=1){    
        ##Combine listed data frames into one dataframe, Labels into one vector
        DataEig <- bind_rows(checkedData)
        LabelsEig <- unlist(checkedLabels)
	rownames(DataEig) <- names(LabelsEig)
	wgRes <- combine.networks(nets=nets, contributions=cont, outPath=saveDir,     
                                  RsquaredCut=RsquaredCut, minModuleSize=20,   
                                  datExpr=DataEig)

    } else {
	DataEig <- DataI
	LabelsEig <- LabelsI
    }

    extractedIDs <- lapply(DataEig, rownames)
    if(length(Reduce(intersect, extractedIDs))!=0)
        stop("Cannot have same row ID in multiple data sets.")
    results[["moduleRes"]] <- wgRes
## Eigengenes:
    pigengene <- compute.pigengene(Data=DataEig, Labels=LabelsEig,
                                   modules=wgRes$modules,
                                   saveFile=combinedPath(saveDir, 'pigengene.RData'),
                                   doPlot='TRUE', verbose=verbose, dOrderByW=dOrderByW)
    results[["pigengene"]] <- pigengene
    ## Multiple conditions?
    if(length(unique(unlist(Labels)))<2){
        message.if("A single condition in Labels. BN and trees are skipped.",
                   verbose=verbose)
        return(results)
    }
    selectedFeatures <- NULL
    ## BN:
    if(bnNum !=0){
        ## Arguments:
        bnPath <- combinedPath(saveDir, 'bn')
        bnArgs$bnPath <- if(!is.null(bnArgs$bnPath)) bnArgs$bnPath else bnPath
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
    c5treeRes <- make.decision.tree(pigengene=pigengene, Data=DataEig,
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
