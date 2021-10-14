consensus <- function(
    candidates, ratio=1/3, bnInputFile=NULL, 
    threshold="Auto", saveFile=NULL, candlist=NULL, bnet.indv=NULL, Data, 
    blacklist, doRemoveBlk=TRUE, verbose=0) 
{
    ## if threshold="Auto", (mean + sd) of the strenghts will be used. 
    message.if("Computing consensus...", verbose=verbose)
    if(is.null(Data))
        stop("Data cannot be NULL")
    scores <- as.numeric(candidates[, "Score"])
    mq <- quantile(scores, (1 - ratio))
    whichcands <- candidates[which(scores >=mq), , drop=FALSE]
    num <- nrow(whichcands)
    top.bnet.indv <- list()
    for (i in 1:num) {
        file1 <- whichcands[i, "File"]
        if (!file.exists(file1)) 
            stop(paste("No file at:", file1, " Bad candidates input! Recompute."))
        bnets <- get(load(file1)) ##bnets
        bnet.indv <- bnets$bnet.indv
        top.bnet.indv[[i]] <- bnet.indv[[as.numeric(whichcands[i, "Index"])]]
        if(doRemoveBlk)
            top.bnet.indv[[i]] <- enforce.blk(bnStrength=top.bnet.indv[[i]], 
                                              blacklist=blacklist)$cleaned
    }
    f1 <- custom.strength(top.bnet.indv, colnames(Data), cpdag=FALSE)
    ##^ Without cpdag=FALSE, edges to Disease may occur. --Habil.
    if(threshold=="Auto"){
        threshold <- mean(f1[, "strength"])+sd(f1[, "strength"])
        threshold <- min(threshold, 1)
        message.if(paste("threshold is automatically determined to be:", threshold), 
                   verbose=verbose)
    }
    c1 <- averaged.network(strength=f1, threshold=threshold)
    c12 <- pdag2dag(c1, ordering=colnames(Data))
    sc <- score(c12, Data, "bde")
    scBIC <- score(c12, Data, "bic")
    message.if(paste("**** The consensus network of", num, "bnet.indvs has a BDe score of ", 
                     sc, " and a BIC score of ", scBIC, "\n"), verbose=verbose)
    consensusRes <- list()
    consensusRes[["threshold"]] <- threshold
    consensusRes[["bnet"]] <- f1
    consensusRes[["BN"]] <- c12
    consensusRes[["PDAG"]] <- c1
    consensusRes[["BDe"]] <- sc
    consensusRes[["BIC"]] <- scBIC
    consensusRes[["Data"]] <- Data
    consensusRes[["num"]] <- num
    consensusRes[["saveFile"]] <- saveFile
    consensusRes[["blacklist"]] <- blacklist
    if (!is.null(saveFile)) {
        save.if(consensusRes, file=saveFile, verbose=verbose)
        ##message.if(paste("consensusRes was saved in:", saveFile))
    }
    return(consensusRes)
}
