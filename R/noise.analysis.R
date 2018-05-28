noise.analysis <- function(
    c5tree, pigengene, Data, testD=NULL, 
    testL=NULL, trainTypes=NULL, saveDir=".", noise, seed=NULL, jump="Auto", verbose=0, 
    repNum=10000, xlab1="Noise (% of expression profile)")
{
    ## Runs add.noise repNum times & takes the average of accuracy for each noise level.
    ## noise: A value in [0, 1]. If not 0, upto this portion of the test expression
    ## will be replaced by Gaussian noise to estimate sensitivity to noise.
    ## jump: The number of entries manuplated in each iteration.
    ##
    if(verbose>0)
        message("Noise analysis...")
    res <- list()
    set.seed(seed)
    dir.create(path=saveDir, recursive=TRUE, showWarnings=FALSE)
    saveFile <- combinedPath(saveDir, "noisy.RData")
    plotFile <- combinedPath(saveDir, "noisy.png")
    accuracyMatrix <- c()
    inds <- 1:repNum
    for(ind in inds){
        if(verbose>1)
            message(ind)
        added <- add.noise(c5tree=c5tree, pigengene=pigengene, 
                           Data=Data, testD=testD, 
                           testL=testL, trainTypes=trainTypes, noise=noise, 
                           seed=NULL, jump=jump, verbose=verbose-1)
        accuracyMatrix <- rbind(accuracyMatrix, added$accuracies)
    }
    rownames(accuracyMatrix) <- inds
    ## Plot:
    x1 <- as.numeric(colnames(accuracyMatrix))
    y1 <- colMeans(accuracyMatrix)
    sd1 <- colSds(accuracyMatrix)
    png(plotFile, height=2*480, width=2*480)
    par(mar=par()$mar+c(1, 1, 0, 0))
    plot(x=x1, y=y1, col='red', ylim=range(c(y1-sd1, y1+sd1)), ylab="Accuracy", xlab=xlab1, 
         cex.axis=2, cex.lab=3, cex=3, pch=20)
    arrows(x0=x1, y0=y1-sd1, x1=x1, y1=y1+sd1, code=3, length=0.1, angle=90, col='grey')
    dev.off()
    ## Output:
    res[['accuracyMatrix']] <- accuracyMatrix
    res[['seed']] <- seed
    res[['jump']] <- jump
    res[['saveDir']] <- saveDir
    res[['saveFile']] <- saveFile
    noisy <- res
    save(noisy, file=saveFile)
    if(verbose>0)
        message(paste("noisy was saved at:", saveFile))
    return(res)
}
