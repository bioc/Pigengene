add.noise <- function(
    c5tree, pigengene, Data, testD=NULL, 
    testL=NULL, trainTypes=NULL, noise, seed=NULL, jump="Auto", verbose=2)
{
    ## Adds Gaussian noise to the testD and reports how the accuracy is affected.
    ## noise: A value in [0, 1]. If not 0, upto this portion of the test expression
    ## will be replaced by Gaussian noise to estimate sensitivity to noise.
    ## jump: The number of entries manuplated in each iteration.
    if(verbose>0)
        message("Adding noise")
    res <- list()
    set.seed(seed)
    accuracies <- c()
    ## Data:
    if(is.null(testD) & is.null(testL)){
        testD <- Data
        if(is.null(trainTypes))
            stop("I cannot compute accuracy without trainTypes and testD!")
        testL <- trainTypes
    }
    N <- length(testD) ## Total number of entries.
    mu <- mean(testD)
    ro <- sd(testD)
    if(verbose>0){
        message(paste("mu=", mu))
        message(paste("rho=", ro))
    }
    if(jump=="Auto")
        jump <- N %/% 100
    noiseNum <- 0 ## Number of noise enteries in testD.
    while(noiseNum/N < noise){ ## Not too much noise, 
        percentage <- 100*(noiseNum/N)
        if(verbose>1)
            message(paste("noiseNum=", noiseNum, ", ", percentage, "%"))
        ## Add noise
        noisyData <- testD
        noisyData[sample(size=noiseNum, x=1:N)] <- rnorm(n=noiseNum, mean=mu, sd=ro)
        ## Projecting:
        p1 <- project.eigen(Data=noisyData, pigengene=pigengene, verbose=0)
        ## Predicting:
        t1 <- table(predict(c5tree, as.data.frame(p1$projected)), testL)
        accuracies[as.character(percentage)] <- sum(diag(t1))/sum(t1)
        noiseNum <- noiseNum+jump
    }
    res[['accuracies']] <- accuracies
    res[['seed']] <- seed
    res[['jump']] <- jump
    return(res)
}
