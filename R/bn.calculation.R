bn.calculation <- function(
    bnInputFile, moduleNum=60, numClust=2, 
    perJob=10, algo='hc', scoring='bde', 
    indivBNs=TRUE, restart=0, pertFrac=0.1, pertFrac0=0, 
    resultPath, ind='', seed=NULL, bnStartFile="None", start=NULL, 
    doRevblk=FALSE, randomInit=TRUE, revBlkRate=0, selected=NULL, 
    Data=NULL, blacklist="Auto", verbose=0, 
    maxSeconds=5*60, doShuffle=TRUE)
{
    ## Feel free to edit, but please do not unnecessarily remove comments! --Amir
    ## bnInputFile is the output of make.bn.input.R, 
    ## bnInputFile contains bnInp, 
    ## a list that has Data (discretesized expression, and optionally Disease, beAML, dummies) blacklist
    ## if indivBNs==T, bn.boot is used and individual networks are saved along with a strength summary, 
    ## otherwise boot.strength is used and only the strength summary is saved. 
    ## If restart > 0, then at each restart, hc's perturb argument will be set
    ##to (pertFrac* (number of nodes)).
    ## ind helps with the future implementation of wrapper/harvester function, use numbers
    ## bnStartFile: contanins the starting BN. Set to "None" to disable.
    ## doRevblk: If TRUE and bnStartFile is "None", the reverse of all edges in the blacklist
    ## will be in the starting (initial) network.
    ## randomInit: Boolean; start from a random network; ignored if either start or bnStartFile are given.
    ## pertFrac0: Is start is provided, its arcs will be perturbated (removed, pruned) by this portion.
    ## maxSeconds: The maximum time allowed for learning the bn if indivBNs=FALSE.
    ## For a network with 30-40 nodes, maxSeconds=5*60 is recommended. 
    ## Beyond this limit, the computation is intrupted and nothing is saved.
    ## If it is not Inf, only 1 core will be used, even if numClust > 1.
    ## blacklist: Set to "Auto" to read from the bnInp list in bnInputFile.
    result <- list()
    files <- c()
    startTime <- Sys.time()
    message.if(paste("Learning", perJob, " BNs using", algo, " started at:"), verbose=verbose)
    message.if(as.character(startTime), verbose=verbose)
    message.if(paste("doRevblk:", doRevblk), verbose=verbose)
    result[["timedOut"]] <- FALSE
    ##
    result[["seed"]] <- seed
    if(!is.null(seed)){
        set.seed(seed)
        message.if(paste("Random seed is set to:", seed), verbose=verbose)
    }
    bnInp <- get(load(bnInputFile)) ## bnInp
    ##^ bnInputFile contains bnInp, 
    ##which has Data (discretesized expression)
    ##and optionally blacklist for Disease, beAML, and dummies.
    if(is.null(Data))
        Data <- bnInp$Data
    if(is.null(Data))
        stop("Data cannot be NULL!")
    if(blacklist=="Auto")
        blacklist <- bnInp$blacklist 
    if(is.null(blacklist))
        warning("blacklist is NULL!")
    if(!bnStartFile=="None"){
        if(!is.null(start))
            stop("Please set either start to NULL, or bnStartFile to 'None'. Too many inputs!")
        load(bnStartFile) ## selected
        start <- selected$BN
    } else {
        if(doRevblk){ ## Add the reverse of blacklist in the start network
            nodes1 <- colnames(Data)
            arcs1 <- blacklist[, c("to", "from")]
            colnames(arcs1) <- c("from", "to")
            start <- empty.graph(nodes=nodes1)
            arcs(start) <- arcs1
        }
    }## Enf if(!bnStartFile=="None").
    
    ## Initial network:
    message.if("Preparing the initial network...", verbose=verbose-1)
    w1 <- "IGNORING randomInit , because start or bnStartFile had been provided !!"
    if(randomInit==TRUE){
        if(!is.null(start) ){
            warning(w1)
        } else { ## no start, 
            nodesOrdered <- colnames(Data)
            if("Disease" %in% nodesOrdered)
                ## Move Disease to the begining. 
                nodesOrdered <- c("Disease", nodesOrdered[-which(nodesOrdered=="Disease")])
            start <- sample.network(nodes=nodesOrdered, blacklist=blacklist, 
                                    revBlkRate=revBlkRate, doShuffle=doShuffle, 
                                    verbose=verbose-1)
        }
    }
    ## 
    arglist <- list()
    ## specifying the initial network (start) works with hc, 
    ##^but might not work with other learning methods
    if(algo %in% c('hc', 'tabu')){ 
        arglist <- list('blacklist'=blacklist, 'score'=scoring, 'start'=start, 'restart'=restart, 
                        'perturb'=round(pertFrac*ncol(Data)))
    }else{
        arglist <- list('blacklist'=blacklist, 'score'=scoring, 'restart'=restart , 
                        'perturb'=round(pertFrac*ncol(Data)))
    }##end algo is hc
    ##
    if(!is.null(start)){
        result[["bdeBefore"]] <- score(start, data=Data, type="bde")
        message.if("BDe score before trainig:", verbose=verbose-1)
        message.if(result[["bdeBefore"]], verbose=verbose-1)
        ## Perturbation:
        arcNum <- nrow(arcs(start))
        pertNum <- round(pertFrac0*arcNum)
        if(pertNum>0){
            arcs(start) <- arcs(start)[-sample(size=pertNum, 1:arcNum, replace=FALSE), ]
            result[["bdeBeforePerturb"]] <- score(start, data=Data, type="bde")
            message.if("BDe score before trainig, perturbed:", verbose=verbose-1)
            message.if(result[["bdeBeforePerturb"]], verbose=verbose-1)
        }
    }
    result[["realStart"]] <- start
    result[["arglist"]] <- arglist
    clustOrNot <- NULL
    if(maxSeconds==Inf){
        clustOrNot <- parallel::makeCluster(numClust)
        ##^ We can change the number of clusters. On stampede numClust can be 16-32 or more
        ## To fix the warning related to makeCluster
        fixArgs <-  c("cluster", "MPIcluster", "PVMcluster", "SOCKcluster")
        assignInNamespace("supported.clusters",fixArgs, "bnlearn")
    }
    ## BN Calculation
    if(! indivBNs ){
        bnet <- boot.strength(data=Data, R=perJob, algorithm=algo, cluster=clustOrNot, 
                              algorithm.args=arglist )
    }else{
        limited.bn.boot <- function(maxSeconds){ ## To disallow long time computation.
            setTimeLimit(cpu=maxSeconds, transient=TRUE)
            bb <- bn.boot(data=Data, R=perJob, algorithm=algo, statistic=arcs, 
                          algorithm.args=arglist, cluster=clustOrNot)
            setTimeLimit(cpu=Inf, elapsed=Inf)
            return(bb)
        }
        tried <- try(bnet.indv <- limited.bn.boot(maxSeconds=maxSeconds), silent=TRUE)
        if(inherits(tried, "try-error")){ ## Timed out.
            m1 <- "Computation of bn.boot() was terminated as its running time exceeded"
            message.if(paste(m1, maxSeconds, "seconds."), verbose=verbose)
            ## Clean up
            if(!is.null(clustOrNot))
                parallel::stopCluster(clustOrNot)
            result[["timedOut"]] <- TRUE
            return(result)
        }
        ## Clean up
        if(!is.null(clustOrNot))
            parallel::stopCluster(clustOrNot)
        ## save the collection of individual networks
        fn <- indFileName(moduleNum=moduleNum, perJob=perJob, ind=ind, typePhrase="bnet.indv")$names
        files["bnet.indv"] <- combinedPath(resultPath, fn)
        bnets <- list("bnet.indv"=bnet.indv)
        ##message.if(paste("Saving bnets in:", files["bnet.indv"]))
        save.if(bnets, file=files["bnet.indv"], verbose=verbose)
        candlist1 <- cbind(files["bnet.indv"], 1:perJob)
        colnames(candlist1) <- c("File", "Index")
        bnets[["scores"]] <- scoreCandidates(candlist=candlist1, Data=Data, 
                                             scoring='bde', saveFile=NULL, verbose=verbose-1)
        ##save.if(bnet.indv, bnet.indv.scores, file=files["bnet.indv"])
        save.if(bnets, file=files["bnet.indv"], verbose=verbose-1)
        s1 <- summary(as.numeric(bnets[["scores"]]$candidates[, "Score"]))
        message.if(paste("Summary of ", perJob, " scores:"), verbose=verbose-1)
        message.if(s1, verbose=verbose-1)
        ## make a summary 
        bnet <- custom.strength(bnet.indv, colnames(Data))
    }##end else (i.e. indivBNs==TRUE)
    ##
    ## End time:
    message.if("Learning BN took:", verbose=verbose)
    timeTaken <- Sys.time()-startTime
    message.if(as.character(timeTaken), verbose=verbose)
    result[["timeTaken"]] <- timeTaken
    ## save the summary
    fn <- indFileName(moduleNum=moduleNum, perJob=perJob, ind=ind, typePhrase="bnet.strength")$names
    files["bnet"] <- combinedPath(resultPath, fn)
    message.if(paste("Saving bnet in:", files["bnet"]), verbose=verbose)
    result[["files"]] <- files
    result[["algo"]] <- algo
    result[["bnet"]] <- bnet
    bnCalculated <- result
    ##save.if(bnet, learnt, file=files["bnet"])
    save.if(bnCalculated, file=files["bnet"], verbose=verbose)
    ##print("bnCalculated was saved in:", file=files["bnet"])
    return(result)
} ## End of function.
