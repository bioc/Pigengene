balance <- function(Data, Labels, amplification=5, verbose=0, naTolerance=0.05){
    ## Balances Data by oversampling based on Labels so that all types
    ##^have roughly the same number of samples.
    message.if(me="Balancing...", verbose=verbose)
    result <- list()
    result[["call"]] <- match.call()
    if(is.null(names(Labels))){
        names(Labels) <- rownames(Data)
    }
    ## QC:
    c1 <- check.pigengene.input(Data=Data, Labels=Labels, na.rm=TRUE, naTolerance=naTolerance)
    Data <- c1$Data
    Labels <- c1$Labels
    Data <- Data[names(Labels), , drop=FALSE]
    condNames <- unique(Labels)
    for(h1 in condNames){
        assign(paste("nCond", h1, sep=""), sum(Labels==h1))
        dataH1 <- Data[which(Labels==h1), , drop=FALSE]
        assign(paste("Data", h1, sep=""), dataH1)
    }
    myDat <- NULL
    brkpts <- 1
    origSampleInds <- NULL ## The indices of rows
    ##^corresponding to the original samples before balancing.
    Reptimes <- c()
    m2 <- "Oversampling to: "
    for(h1 in condNames){
        ##onm <- paste("cond", h1, "RepTimes", sep='')
        oncond <- get(paste("nCond", h1, sep=''))
        rptim <- round(amplification * (nrow(Data)/oncond))
        if(length(unique(table(Labels)))==1) ## All sampls have the same size,
            rptim <- 1 ## Do not oversample.
        origSampleInds <- c (origSampleInds, brkpts:((brkpts+oncond)-1))
        Reptimes[h1] <- rptim
        m2 <- paste(m2, rptim*oncond, 'of type', h1, ", ")
        repeated <- repeat.data(Data=get(paste("Data", h1, sep="")), times=rptim)$repeated
        myDat <- rbind(myDat, repeated)
        brkpts <- (1+nrow(myDat))
    }
    message.if(me=m2, verbose=verbose)
    result[["origSampleInds"]] <- origSampleInds
    result[['Reptimes']] <- Reptimes
    result[['balanced']] <- myDat
    return(result)
}
