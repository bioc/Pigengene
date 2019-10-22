calculate.beta <- function(saveFile=NULL, RsquaredCut=0.8, Data,  doThreads=FALSE, verbose=0){
    ## RsquaredCut: The default in WGCNA is 0.85 but it may be too much for WGCNA.
    message.if(me="Calculating beta...", verbose=verbose)
    result <- list()
    result[["call"]] <- match.call()
    if(!is.null(saveFile)){
        outFile <- paste(saveFile,".out",sep="")
        fOut <- file(outFile)
        result[["outFile"]] <- outFile
    } else {
        fOut <- file()
    }
    ## Following command tells R to not to convert char variables to factors.
    if(is.null(Data))
        stop("Data cannot be NULL!")
    options(stringsAsFactors=FALSE)
    if(doThreads)
        WGCNA::allowWGCNAThreads()
    powers <- c(c(1:14), seq(from=16, to=20, by=2))    
    sink(file=fOut) ## silence upcoming output 
    ## Condition:
    sft <- pickSoftThreshold(Data, powerVector=powers, 
                             verbose=verbose-1, RsquaredCut=RsquaredCut)
    sink()
    close(fOut)
    ##
    result[['sft']] <- sft
    result[['power']] <- sft$powerEstimate
    message.if(me=paste("beta:", sft$powerEstimate), verbose=verbose)
    result[['powers']] <- powers
    result[['RsquaredCut']] <- RsquaredCut
    calculateBetaRes <- result
    save.if(calculateBetaRes, file=saveFile)
    return(calculateBetaRes)
}
