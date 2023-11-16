scoreCandidates <- function(candlist, scoring="bde", 
                            saveFile=NULL, Data, verbose=0){
    ##
    result <- list()
    message.if("Scoring candidates...", verbose=verbose)
    bestBN <- empty.graph(nodes=colnames(Data))
    failedFiles <- vector(mode="character")
    noAccessFiles <- vector(mode="character")
    files <- unique(candlist[, 1])
    candlist <- candlist[, c("File", "Index")]
    fa <- file.access(files, 0)
    message.if(paste("Number of expected files:", nrow(candlist)), verbose=verbose-1)
    if (sum(fa==-1) > 0) {
        message.if("---- cannot access the following files, files do not exist:", verbose=verbose-1)
        message.if(files[fa==-1], verbose=verbose-1)
        failedFiles <- as.character(files[fa==-1])
        files <- files[fa > -1]
    }
    fa <- file.access(files, 4)
    if (sum(fa==-1) > 0) {
        message.if(">>>> cannot access the following files, permission denied:", verbose=verbose-1)
        message.if(files[fa==-1], verbose=verbose-1)
        noAccessFiles <- as.character(files[fa==-1])
        files <- files[fa > -1]
    }
    candlist <- candlist[(candlist[, 1] %in% files), , drop=FALSE]
    message.if(paste("Number of available files:", nrow(candlist)), verbose=verbose-1)
    scs <- vector(mode="numeric", length=nrow(candlist))
    candlist <- cbind(candlist, scs)
    colnames(candlist)[ncol(candlist)] <- "Score"
    for (i in 1:nrow(candlist)){
        ##load.if(candlist[i, 1], verbose=verbose>2) ## bnets
        candy <- getBN(candlist[i, 1:2], nodelist=colnames(Data))
        rep <- candlist[i, 2]
        message.if(paste("Now scoring: network number ", rep, " from file ", 
                         candlist[i, 1]), verbose=verbose-2)
        candlist[i, 3] <- bnlearn::score(candy, data=Data, type=scoring)
        ## try_with_time_limit( bnlearn::score(candy, data=Data, type=scoring), 4)
    }
    message.if(paste("**** Best Candidate , out of ", nrow(candlist), 
                     " considered candidates is :", "\n"), verbose=verbose-1)
    message.if(candlist[which.max(candlist[, 3]), ], verbose=verbose-1)
    bestBN <- getBN(candlist[which.max(candlist[, 3]), 1:2], colnames(Data), 
                    verbose=verbose-1)
    result[["candidates"]] <- candlist
    result[["failedFiles"]] <- failedFiles
    result[["noAccessFiles"]] <- noAccessFiles
    result[["bestBN"]] <- bestBN
    scoreCandRes <- result
    if (!is.null(saveFile)){
        result[["saveFile"]] <- saveFile
        save.if(scoreCandRes, file=saveFile, verbose=verbose)
    }
    return(scoreCandRes)
}
