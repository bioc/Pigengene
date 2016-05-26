compact.tree <- function(
    c5Tree, pigengene, Data=pigengene$Data, Labels=pigengene$Labels,
    testD=NULL, testL=NULL, saveDir=".", verbose=0)
{
    ##stopifnot((is.null(testD) & is.null(testL)) |
    ##!((is.null(testD) & is.null(testL))), "test data?")
    compRes <- list()
    message.if("Compacting the tree...", verbose=verbose)
    compRes[["call"]] <- match.call()
    dir.create(path=saveDir, recursive=TRUE, showWarnings=FALSE)
    txtFile <- combinedPath(dir=saveDir, fn="compact.txt")
    capture.output(print("A report on tree compacting."), file=txtFile, append=FALSE)
    inTxt <- function(what){
        capture.output(print(what), file=txtFile, append=TRUE)
    }
    if(nrow(Data)==0)
        stop("Data has no rows!")
    if(is.null(testD) & is.null(testL)){
        testD <- Data
        testL <- Labels
    }
    feats <- get.used.features(c5Tree)
    featsN <- as.numeric(gsub(feats, pattern="ME", replacement=""))
    inds <- which(pigengene$orderedModules %in% featsN)
    genes <- names(pigengene$orderedModules[inds])
    inTxt("Modules and number of genes before compacting:")
    inTxt(table(pigengene$orderedModules[genes]))
    modules <- pigengene$orderedModules[colnames(Data)]
    queue <- make.membership.queue(feats, pigengene, modules=modules)
    tmpEig <- pigengene
    qupos3 <- 0
    ERR3 <- 5
    ERR3Ms <- c()
    ERR3s <- c()
    ##colnames(tmpEig$eigengenes) <- colnames(BCR)
    testD <- as.matrix(testD)
    p1 <- project.eigen(Data=testD, pigengene=tmpEig)
    compRes[["predTrain"]] <- predict(c5Tree, as.data.frame(p1$projected))
    startM <- table(compRes[["predTrain"]], testL)
    startErr <- startM[2]+ startM[3]
    ## while(ERR3< (0.4*nrow(testD)) &qupos3< (length(queue)-10) ){
    hTresh <- startErr+(0.1*nrow(testD))
    m1 <- paste(unique(testL), "-err", collapse=" ", sep="")
    message.if("Affect on the test dataset:", verbose=verbose-2)
    message.if(paste("Removed Left", m1, "Total-err"), verbose=verbose-2)
    if(-1 %in% tmpEig$orderedModules)
        stop("Compacting a tree having module -1 is not implemented!") 
    while(ERR3< (hTresh) &qupos3< (length(queue)-10) ){
        modulesTmp <- tmpEig$orderedModules
        curmod <- modulesTmp[names(queue[qupos3])]
        ##^ current module for this position.
        ##
        if(qupos3>0) ## We want to plot without compacting too.
            if((table(modulesTmp)[match(curmod, names(table(modulesTmp)))])<3){
                qupos3 <- qupos3+1 ## the next position.
                next
            }
        tmpEig$orderedModules[names(queue[qupos3])] <- -1
        ## On train:
        P2c3 <- project.eigen(Data=as.matrix(Data), pigengene=tmpEig, ignoreModules=-1)
        ert3M <- table(predict(c5Tree, as.data.frame(P2c3$projected)), Labels)
        ERR3M <- ert3M[2]+ert3M[3]
        ERR3Ms <- rbind(ERR3Ms, c( qupos3, ert3M[2], ert3M[3], ERR3M))
        ## On test  
        BCR2.3 <- project.eigen(Data=testD, pigengene=tmpEig, ignoreModules=-1)
        ert3 <- table(predict(c5Tree, as.data.frame(BCR2.3$projected)), testL)
        ERR3 <- ert3[2]+ert3[3]
        ERR3s <- rbind(ERR3s, c( qupos3, ert3[2], ert3[3], ERR3))
        left <- length(genes)-qupos3
        message.if(paste(qupos3, left, ert3[2], ert3[3], ERR3, collapse=" "), 
                   verbose=verbose-2)
        qupos3 <- qupos3+1 ## the next position.
    }
    ## Genes:
    inds2 <- which(tmpEig$orderedModules %in% featsN)
    genesCompacted <- names(tmpEig$orderedModules[inds2])
    inTxt("Number of genes after compacting:")
    inTxt(table(tmpEig$orderedModules[genesCompacted]))
    ## confusion matrix:
    compRes[["predTrainCompact"]] <- predict(c5Tree, as.data.frame(P2c3$projected))
    inTxt("The confusion matrix from the full tree on train (after projection):")
    ## Replace table with caret::confusionMatrix for reporting performance.
    inTxt(table(compRes[["predTrain"]], Labels))
    inTxt("The confusion matrix from the compacted tree on train:")
    inTxt(table(compRes[["predTrainCompact"]], Labels))
    
    ## Plot error behaviour in response to compacting.
    xlab1 <- 'Number of removed genes'
    ylab1 <- 'Number of misclassified cases'
    mLegend <- c("Total", unique(testL))
    ## Train:
    pngfn <- combinedPath(dir=saveDir, fn=paste('Compacting_test.png', sep=''))
    best1 <- max(which(ERR3s[, 4]==min(ERR3s[, 4])))
    m2 <- paste("Error change on the test dataset \n best on test after removing ", best1)
    png(pngfn)
    plot(ERR3s[, 4], main=m2, ylab=ylab1, xlab=xlab1, ylim=range(ERR3s[, 2:4]))  
    points(ERR3s[, 3], pch=8, col='blue')
    points(ERR3s[, 2], pch=8, col='red')
    abline(v=max(which(ERR3s[, 4]==min(ERR3s[, 4]))), col='pink')
    legend("topleft", col=c("black", "blue", "red"), pch=c(1, 8, 8), mLegend)
    dev.off()
    ## Test:
    best2 <- max(which(ERR3Ms[, 4]==min(ERR3Ms[, 4])))
    m3 <- paste("Error change on the train dataset \n best after removing", best2)
    pngfn <- combinedPath(dir=saveDir, fn=paste('Compacting_train.png', sep=''))
    png(pngfn)
    plot(ERR3Ms[, 4], main=m3, ylab=ylab1, xlab=xlab1, col='orange', ylim=range(ERR3Ms[, 2:4]))
    points(ERR3s[, 3], pch=8, col='blue')
    points(ERR3s[, 2], pch=8, col='red')
    abline(v=max(which(ERR3s[, 4]==min(ERR3s[, 4]))), col='pink')
    legend("topleft", col=c("orange", "blue", "red"), pch=c(1, 8, 8), mLegend)
    dev.off()
    message.if(paste("Changes to accuracy were plotted in", saveDir), verbose=verbose)
    m3 <- "Confusion matrices and other details on compacting were reported in"
    message.if(paste(m3, txtFile), verbose=verbose)
    ##
    trainErrors <- ERR3Ms[,2:4]
    colnames(trainErrors) <- c(unique(testL), "Total")
    rownames(trainErrors) <- ERR3Ms[,1]
    testErrors <- ERR3s[,2:4]
    colnames(testErrors) <- colnames(ERR3Ms)
    rownames(testErrors) <- ERR3s[,1]
    ## Output:
    compRes[['genes']] <- genes 
    compRes[['genesCompacted']] <- genesCompacted 
    compRes[['trainErrors']] <- trainErrors 
    compRes[['testErrors']] <- testErrors
    compRes[['queue']] <-  queue
    compRes[['startErr']] <- startErr
    compRes[['pos']] <- qupos3-1
    compRes[['txtFile']] <- txtFile
    invisible(compRes)
}
