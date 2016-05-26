getALLRuns <- function (dir, moduleNum, perJob, inds=1:100){
    ## returns a comprehensive candidate list, containing every replicate from every Run
    ## use this in combination with scoreCandidates to examine the development of scores.
    ## dir: location of BNs e.g. '/work/03270/zare/proj/genetwork/result/AML/bn/2015-03-02-bnCalc/'
    ## suffix.ind is a vector: set of indices used in the submitted jobs
    ## x below will contain the hypothetical list of results
    ## I.e. what one would expect based on submitted/planned jobs
    result <- list()
    indFileNames <- indFileName(dir=dir, moduleNum=moduleNum, 
                                perJob=perJob, ind=inds)
    x <- combinedPath(dir, indFileNames$names)
    perJob <- indFileNames$perJob
    ca <- cbind(rep(x, each=perJob), rep(c(1:perJob), 
                           length(x)))
    colnames(ca) <- c("File", "Index")
    result[["perJob"]] <- perJob
    result[["candidates"]] <- ca
    allRurnsRes <- result
    return(allRurnsRes)
}
