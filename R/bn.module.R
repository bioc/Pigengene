bn.module <- function(
    resultPath, moduleNum, use.Hartemink=TRUE, bnNum, perJob, 
    restart=0, use.Effect, use.Disease, dummies=NULL, seed=NULL, 
    scoring, pertFrac, bnStartFile="None", 
    doLocal="Auto", saveToPath, selectedFeatures=NULL, bnCalculationJob=NULL, 
    timeJob=NULL, algo="hc", verbose=0, 
    partition=NULL, selected=NULL, Data=NULL, Labels=NULL, 
    pvalGenes=NULL, maxSeconds=5*60, moduleFile="Auto", doShuffle=TRUE, naTolerance=0.05){

    ##
    result <- list()
    message.if(paste("Start learning", bnNum, "BNs for module", moduleNum, "using", algo), 
               verbose=verbose)
    if(moduleFile=="Auto")
        moduleFile <- get.module.data.file(resultPath=resultPath, 
                                           moduleNum=moduleNum, use.Hartemink=use.Hartemink, 
                                           partition=partition)$moduleFile
    message.if("moduleFile:", verbose=verbose)
    message.if(moduleFile, verbose=verbose)
    cluster <- which.cluster()
    queue <- cluster$queue
    if (is.null(timeJob)) 
        timeJob <- cluster$timeJob
    numClust <- cluster$numClust
    maxJobs <- cluster$maxJobs
    if (doLocal=="Auto") 
        doLocal <- !cluster$onCluster
    if (!bnStartFile=="None") {
        load(bnStartFile)
        pertFrac <- selected$pertFrac
        result[["pertFrac"]] <- pertFrac
    }
    made <- make.bn.input(moduleNum=moduleNum, use.Hartemink=use.Hartemink, 
                          use.Effect=use.Effect, use.Disease=use.Disease, breakNum=3, 
                          ibreaks=20, dummies=dummies, saveFile=moduleFile, 
                          selectedFeatures=selectedFeatures, verbose=verbose-1, 
                          Data=Data, Labels=Labels , pvalGenes=pvalGenes , naTolerance=naTolerance)
    bnS <- bn.calculationS(moduleNum=moduleNum, bnPath=saveToPath, verbose=verbose-2, 
                           totalNum=bnNum, bnCalculationJob=bnCalculationJob, 
                           pertFrac=pertFrac, moduleFile=moduleFile, doLocal=doLocal, 
                           queue=queue, scoring=scoring, perJob=perJob, restart=restart, 
                           timeJob=timeJob, numClust=numClust, maxJobs=maxJobs, 
                           seed=seed, bnStartFile=bnStartFile, inds="Auto", 
                           algo=algo, maxSeconds=maxSeconds, doShuffle=doShuffle)
    result[["moduleFile"]] <- moduleFile
    ##^ The output of make.bn.input() containing data and blacklist arcs.
    result[["bnS"]] <- bnS
    result[["algo"]] <- algo
    result[["seed"]] <- seed
    bnModuleRes <- result 
    return(bnModuleRes)
}
