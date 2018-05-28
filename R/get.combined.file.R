get.combined.file <- function (ind, moduleNum, perJob, resultPath) 
{
    res <- list()
    mimx <- paste(min(ind), max(ind), sep="-")
    resfn <- indFileName(moduleNum=moduleNum, perJob=perJob, 
                         ind=mimx, typePhrase="bnet.combined.strength")$names
    res$combinedFile <- combinedPath(resultPath, resfn)
    return(res)
}
