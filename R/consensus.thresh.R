consensus.thresh <- function(
    candidates, ratio=1/3, Data, blacklist, 
    threshVector=c(0.4, 0.5, 0.6), consensusFile=NULL, verbose=0)

{
    result <- list()
    files <- c()
    consensusRes <- list()
    m1 <- "Computing consensus network for threshold in "
    message.if(paste(m1, threshVector), verbose=verbose)
    for (tI in threshVector) {
        tChar <- as.character(tI) ## tIChar
        files[tChar] <- get.consensus.files(
            consensusFile=consensusFile, threshVector=tI, verbose=verbose-1)
        c1 <- consensus(candidates=candidates, blacklist=blacklist, Data=Data, 
                        ratio=ratio, threshold=tI, saveFile=files[tChar], 
                        verbose=verbose-1)
        consensusRes[[tChar]] <- c1
    }
    result[["files"]] <- files
    result[["consensusRes"]] <- consensusRes
    result[["threshVector"]] <- threshVector
    consensusThreshRes <- result
    return(consensusThreshRes)
}
