learn.bn <- function(
    pigengene=NULL, Data=NULL, Labels=NULL, bnPath="bn", bnNum=100, consensusRatio=1/3, 
    consensusThresh="Auto", doME0=FALSE, selectedFeatures=NULL, trainingCases="All", 
    algo="hc", scoring="bde", restart=0, pertFrac=0.1, doShuffle=TRUE, 
    use.Hartemink=TRUE, bnStartFile="None", use.Disease=TRUE, 
    use.Effect=FALSE, dummies=NULL, tasks="All", 
    onCluster=!(which.cluster()$cluster=="local"), 
    inds=1:ceiling(bnNum/perJob), perJob=2, 
    maxSeconds=5*60, timeJob="00:10:00", bnCalculationJob=NULL, 
    seed=NULL, verbose=0)
{
    ## Sets the parameters and learns BN for a module. --Habil.
    ## trainingCases: are the patients that are used for training the model.
    ##^ Set to "All" to use all.
    ## pigengene: The result of compute.pigengene().
    ##^ Used to be eigengeneHeat before Pigengene.
    ## bnCalculationJob: The script that sbatch() needs to submit the jobs to the cluster.
    ## consensusThresh: A numeric vector with values in [0, 1], or set it to "Auto"
    ##^ to automatically choose the threshould as (mean+sd) of strengths. 
    ##^ While making the consensus network, any connection with strength less than
    ##^ this threshoulds will be ignored. If more than 1 value is provided, 
    ## a consensus network will be build for each value. 
    ## maxSeconds: If a training job takes more than this many seconds, 
    ##^it will be terminated. Not effective on cluster.
    ## seed: The seed that is used to produce the random initial network
    ##^and also during the learning procedure.
    ##^ Set to NA to disable.
    ## selectedFeatures: If not NULL, the network will be learnt only on these
    ##^ variables (modules).
    startTime <- Sys.time()
    message.if(paste("learn.bn() with bnNum=", bnNum, "started at:"), verbose=verbose)
    message.if(as.character(startTime), verbose=verbose)
    res <- list() ## result.
    moduleNum <- "E"
    dir.create(path=bnPath, recursive=TRUE, showWarnings=FALSE)
    res[["consensusThresh"]] <- consensusThresh
    res[["call"]] <- match.call()
    combiningThresh <- consensusThresh[1]
    ##^ A numeric in [0, 1] used to make the consensus network from
    ##^ALL the individual networks.
    ## What to do:
    tasksFull <- c("learn", "harvest", "consensus", "toLocal", "graph")
    if("All" %in% tasks)
        tasks <- c("learn", "harvest", "consensus", "graph")
    indvPath <- combinedPath(bnPath, 'indv')
    moduleFile <- get.module.data.file(resultPath=bnPath, 
                                       moduleNum=moduleNum, use.Hartemink=use.Hartemink, 
                                       partition=NULL)$moduleFile
    scoreFile <- combinedPath(bnPath, "allTestScores.RData")
    consensusFile <- combinedPath(bnPath, "consensus.RData")
    res[["indvPath"]] <- indvPath
    res[["moduleFile"]] <- moduleFile
    res[["scoreFile"]] <- scoreFile
    res[["consensusFile"]] <- consensusFile 
    res[["seed"]] <- seed 
    ##
    stop1 <- "Provide an appropriate pigengene, or Data and Labels as input."
    if("learn" %in% tasks){
        if(moduleNum=="E"){ ## eigengenes ##Always TRUE
            ## data:
            if(!is.null(pigengene) & is.null(Data) & is.null(Labels)){
                Data <- pigengene$eigengenes
                Labels <- pigengene$annotation[, 1]
                names(Labels) <- rownames(pigengene$annotation)
            } else {
                if(is.null(Data)|is.null(Labels))
                    stop(stop1)
            }
            c1 <- check.pigengene.input(Data=Data, Labels=Labels, na.rm=TRUE)
            Data <- c1$Data
            Labels <- c1$Labels
            ## Module 0:
            if(!doME0)
                Data <- Data[, colnames(Data)!="ME0"]
            if(trainingCases!="All"){
                Data <- Data[trainingCases, ]
                Labels <- Labels[trainingCases, 1,drop=FALSE]
            }
            selectedFeatures <- colnames(Data)
        }##End if(moduleNum=="E").
        ## Learning:
        dir.create(path=indvPath, recursive=TRUE, showWarnings=FALSE)
        res[["bnModuleRes"]] <- bn.module(
            resultPath=bnPath, saveToPath=indvPath, 
            moduleNum=moduleNum, use.Hartemink=use.Hartemink, 
            bnNum=bnNum, perJob=perJob, 
            use.Effect=use.Effect, use.Disease=use.Disease, 
            dummies=dummies, pertFrac=pertFrac, scoring=scoring, 
            Data=Data, Labels=Labels, 
            bnStartFile=bnStartFile, 
            selectedFeatures=selectedFeatures, 
            bnCalculationJob=bnCalculationJob, timeJob=timeJob, 
            algo=algo, seed=seed, verbose=verbose-1, 
            partition=NULL, maxSeconds=maxSeconds, 
            moduleFile=moduleFile, doShuffle=doShuffle)
        if(onCluster){
            Sys.sleep(180)
        }
    }##End if("learn" %in% tasks).

    bnInp <- get(load(moduleFile)) ## bnInp
    Data <- bnInp$Data
    blacklist <- bnInp$blacklist
    if("harvest" %in% tasks){
        message.if("Harvesting...", verbose=verbose-1)
        res[["runs"]] <- getALLRuns(dir=indvPath, moduleNum=moduleNum, perJob="Auto", 
                                    inds=inds)
        res[["scores"]] <- scoreCandidates(candlist=res$runs$candidates, 
                                           Data=Data, verbose=verbose-1, 
                                           scoring=scoring, saveFile=scoreFile)
        ##^ will score each individual network, note which jobs did not
        ##^finish and which files where not accessible
    }
    
    scoreCandRes <- get(load(scoreFile)) ## scoreCandRes
    candidates <- scoreCandRes$candidates
    if("consensus" %in% tasks){
        ctr <- consensus.thresh(candidates=candidates, ratio=consensusRatio, 
                                Data=Data, blacklist=blacklist, verbose=verbose-1, 
                                threshVector=consensusThresh, 
                                consensusFile=consensusFile)
        res[["consensusThreshRes"]] <- ctr
        consensus1 <- ctr$consensusRes[[as.character(ctr$threshVector[1])]]
        res[["consensus1"]] <- consensus1
    }
    
    if("toLocal" %in% tasks){ ## Copy from cluster to local and draw.scores.
        stop("toLocal task is not implemented yet!")
        ## First copy it to the local machine. E.g.:
        S <- NULL ## This task is under development and not ready for user.
        fromFiles <- gsub(S$savePath, pattern=S$localPath, replacement=S$maverickPath)
        fromFiles <- paste(paste(fromFiles, 
                                 c("/allTestScores.RData", "/*combined*", 
                                   "/consensus/consensus*.RData"), sep=""), collapse=" ")
        command <- paste("scp ", S$clusterUsername, "@", S$clusterUrl, ':"', fromFiles, '" ', 
                         S$savePath, sep="")
        message.if(paste("Enter password for:", command), verbose=1)
        system(command)
        commandMv <- paste("mv ", S$savePath, "/consensus*.RData ", S$consensusPath, sep="")
        system(commandMv) ## Move to consensus folder.
    }

    if("graph" %in% tasks){
        res[["scorePlot"]] <- draw.scores(
            candidates=candidates, verbose=verbose-1, 
            savePath=dirname(scoreFile))
        res[["graphs"]] <- draw.bnS(
            consensusFile=consensusFile, verbose=verbose-1, 
            consensusThresh=consensusThresh)
    }
    ##
    message.if("learn.bn() took:", verbose=verbose)
    timeTaken <- Sys.time()-startTime
    message.if(format(timeTaken), verbose=verbose)
    res[["timeTaken"]] <- timeTaken
    ##
    return(res)
}##End learn.bn <- function.
