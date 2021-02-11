bn.calculationS <- function(
    moduleNum, bnPath, perJob=15, totalNum, doTalk=verbose>0, 
    bnCalculationJob=NULL, inds="Auto", moduleFile, restart, pertFrac=0.1, 
    doLocal=FALSE, queue="normal", scoring="bde", timeJob="47:00:00", 
    numClust, maxJobs, seed=NULL, bnStartFile="None", algo="hc", selected=NULL, 
    maxSeconds=5*60, doShuffle=TRUE, verbose=0)
{
    result <- list()
    if (inds=="Auto") {
        inds <- 1:ceiling(totalNum/perJob)
    }
    files <- list()
    sbatchS <- list()
    set.seed(seed)
    ## We need a different seed for each calculation.
    seeds <- as.integer(runif(n=length(inds), min=1, max=10^9))
    ##^ A seed should be an integer <2*10^9.
    for (ind in inds) {
        seedInd <- seeds[ind]
        verboseInd <- verbose-1
        args <- c(paste("--bnInputFile", moduleFile), 
                  paste("--moduleNum", moduleNum), 
                  paste("--numClust", numClust), 
                  paste("--perJob", perJob), 
                  paste("--algo", algo), 
                  paste("--scoring", scoring), 
                  paste("--indivBNs", TRUE), 
                  paste("--restart", restart), 
                  paste("--pertFrac", pertFrac), 
                  paste("--resultPath", bnPath), 
                  paste("--ind", ind), 
                  paste("--bnStartFile", bnStartFile), 
                  paste("--seed", seedInd), 
                  paste("--doShuffle", doShuffle), 
                  paste("--verbose", verboseInd)
                  )
        if (scoring !="bde") 
            warning("--scoring is NOT 'bde', ## should be 'bde'")
            if (doTalk) {
                message("The arguments are set as:")
                message(paste(args, collapse=" "))
            }
        if (doLocal) {
            m1 <- matrix(unlist(strsplit(args, split=" ")), byrow=FALSE, nrow=2)
            ##^ Does not work if moduleFile contains space!
            colnames(m1) <- gsub(m1[1, ], pattern="--", replacement="")
            opt <- as.list(m1[2, ])
            bns <- bn.calculation(bnInputFile=opt$bnInputFile, 
                                  moduleNum=opt$moduleNum, 
                                  numClust=as.numeric(opt$numClust), 
                                  perJob=as.numeric(opt$perJob), algo=opt$algo, 
                                  scoring=opt$scoring, indivBNs=as.logical(opt$indivBNs), 
                                  restart=as.numeric(opt$restart), 
                                  pertFrac=as.numeric(opt$pertFrac), 
                                  resultPath=opt$resultPath, ind=opt$ind, 
                                  seed=seedInd, bnStartFile=bnStartFile, 
                                  maxSeconds=maxSeconds, 
                                  doShuffle=doShuffle, verbose=as.numeric(opt$verbose))
            files[[as.character(ind)]] <- bns[["files"]]
        } else { ## on cluster
            if(is.null(bnCalculationJob))
                stop("A bnCalculationJob script should be provided to submit the jobs to the cluster!")
            s1 <- sbatch(sourceFile=bnCalculationJob, doSubmit=TRUE, 
                         doDeleteTempFile=FALSE, 
                         jobName=paste("bn", ind, "m", moduleNum, sep=""), doTalk=TRUE, 
                         queue=queue, arguments=args, time1=timeJob, 
                         maxJobs=maxJobs)
            sbatchS[[as.character(ind)]] <- s1
        }
    }
    result[["files"]] <- files
    result[["sbatchS"]] <- sbatchS
    result[["algo"]] <- algo
    result[["seed"]] <- seed
    bnCalculationSRes <- result
    return(bnCalculationSRes)
}
