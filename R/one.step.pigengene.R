one.step.pigengene <- function(Data, saveDir="Pigengene",
                               Labels, testD=NULL, testLabels=NULL, doBalance=TRUE, RsquaredCut=0.8,
                               costRatio=1, toCompact=FALSE, bnNum=0, bnArgs=NULL, useMod0=FALSE,
                               mit="All", ## unique(Labels)[1],
                               verbose=0, doHeat=TRUE, seed=NULL, dOrderByW=TRUE, naTolerance=0.05, doNetOnly=FALSE){
    ## costRatio: Implemented only for 2 classes.
    ##^Determines how severe it is to misclassify a sample accross types.
    ##^E.g., if costRatio=2, misclassification of a sample of the 1st type is
    ##considered twice worse than misclassification of a sample of the 2nd type.

    results <- list()
    results[["call"]] <- match.call()

    dataNum <- 1
    if(inherits(Data, "list")){
        dataNum <- length(Data)
    } else {
        m1 <- paste("Pigengene started analizing", nrow(Data),
                    "samples using", ncol(Data), "genes...")
        message.if(me=m1, verbose=verbose)
        if(verbose>1){
            print(table(Labels))
        }
    }

    ## saveDir:
    if(length(grep(saveDir, pattern=" ")>0))
        stop("saveDir cannot have space!")
    dir.create(path=saveDir, recursive=TRUE, showWarnings=FALSE)

    if(is.null(testD) & !is.null(testLabels)){
        warning("testLabels is ignored because testD is NULL.")
        testLabels <- NULL
    }
    if(!is.null(testD)){
        if(is.null(testLabels))
            stop("testLabels is NULL while testD is not NULL!")
        ct <- check.pigengene.input(Data=testD, Labels=testLabels, na.rm=TRUE, naTolerance=naTolerance)
        testD <- ct$Data
        testLabels <- ct$Labels
    }

    nets <- list()
    cont <- c()
    checkeData <- list()    
    checkedLabels <- list()    
    for(ind in 1:dataNum){
        if(dataNum==1){
	    DataI  <- Data
            LabelsI <- Labels
        } else {
            DataI <- Data[[ind]]
            LabelsI <- Labels[[ind]]
            ## contribution from the given data set DataI
	    cont[ind] <- nrow(DataI)
        }
        ## QC: 
	c1 <- check.pigengene.input(Data=DataI, Labels=LabelsI, na.rm=TRUE, naTolerance=naTolerance)
        DataI <- c1$Data
        LabelsI <- c1$Labels
	checkeData[[ind]] <- as.data.frame(DataI)
	checkedLabels[[ind]] <- LabelsI
        ## Data for WGCNA:
        wData <- switch(mit,
                        "All"=DataI,
                        DataI[which(LabelsI %in% mit), ])
        ##stop('mit must be equal to "cond1", "cond2", or "Both"!'))
        ## WGCNA:
        if(doBalance)
            wData <- balance(Data=wData, Labels=LabelsI, verbose=verbose-1, naTolerance=naTolerance)$balanced
        message.if("Computing correlation...", verbose=verbose-1)
        nets[[ind]] <- abs(stats::cor(wData))
        if(dataNum==1){
            print("dataNum==1")
            calculateBetaRes <- calculate.beta(saveFile=NULL, RsquaredCut=RsquaredCut,
                                               Data=wData, verbose=verbose-1)
            results[["betaRes"]] <- calculateBetaRes
            if(is.na(calculateBetaRes[["power"]]))
                stop("Consider a lower value for RsquaredCut, power is NA!")
            betaI <- calculateBetaRes[["power"]]
            wgRes <- wgcna.one.step(Data=wData, seed=seed,
                                    power=betaI,
                                    saveDir=saveDir, verbose=verbose-1)
            results[["moduleRes"]] <- wgRes
            modules <- wgRes$modules
            ##browser()
            if(doNetOnly){
                warning("Identify modules took time, but is not needed when doNetOnly==TRUE")
            }
            results[["netMatrix"]] <- nets[[ind]]
        }
        rm(wData)
    }
    
    ## Now use combine.networks()
    if(dataNum!=1){
        print("dataNum!=1")
        ##Combine listed data frames into one dataframe, Labels into one vector
        message.if("Binding data...",  verbose=verbose-2)
        DataEig <- as.matrix(dplyr::bind_rows(checkeData))
        LabelsEig <- unlist(checkedLabels)
	## check for duplicates
	extractedIDs <- rownames(DataEig)
	if(any(duplicated(extractedIDs)))
	    stop("Cannot have same row ID in multiple data sets.")
	rownames(DataEig) <- names(LabelsEig)
	combined <- combine.networks(nets=nets, contributions=cont, outPath=saveDir,     
                                  RsquaredCut=RsquaredCut, minModuleSize=20,   
                                  datExpr=DataEig, verbose=verbose-1, doReturNetworks=doNetOnly,
                                  doIdentifyModule=!doNetOnly)
        results[["netMatrix"]] <- combined$netMatrix
        results[["combined"]] <- combined
        if(!doNetOnly){
            modules <- combined$modules
        }
         
    } else { ##dataNum=1
	DataEig <- DataI
	LabelsEig <- LabelsI
    }
    ##browser()
    if(doNetOnly){
        return(results)
    }
    
    ## Eigengenes:
    pigengene <- compute.pigengene(Data=DataEig, Labels=LabelsEig,
                                   modules=modules,
                                   saveFile=combinedPath(saveDir, 'pigengene.RData'),
                                   doPlot='TRUE', verbose=verbose, dOrderByW=dOrderByW,
                                   naTolerance=naTolerance)

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
        bnArgs <- c(bnArgs, list(pigengene=pigengene, bnNum=bnNum, verbose=verbose-1, seed=seed,
                                 naTolerance=naTolerance))
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
        if(length(selectedFeatures)<2){
            warning("The condition variable has <2 children. BN results will be ignored.")
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
                                    toCompact=toCompact, naTolerance=naTolerance)
    ##
    results[["selectedFeatures"]] <- selectedFeatures
    results[["c5treeRes"]] <- c5treeRes
    return(results)
}
