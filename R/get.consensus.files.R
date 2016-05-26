get.consensus.files <- function (consensusFile, threshVector, verbose=0){
    message.if("Getting consensus files...", verbose=verbose)
    message.if(paste("consensusFile:", consensusFile), verbose=verbose-1)
    if (length(threshVector)==0) {
        stop("No threshold given!")
    }
    consensusFiles <- c()
    for (threshI in threshVector) {
        tempConsFile <- gsub(consensusFile, pattern="\\.RData", replacement=paste("-T", 
                                                                    threshI, ".RData", sep=""))
        consensusFiles <- c(consensusFiles, tempConsFile)
    }
    return(consensusFiles)
}
