combinedStrength <- function(dir, bnInputFile=NULL, moduleNum, repsPerRun, inds,
                             typePhrase="bnet.strength", resdir=".", doSave=TRUE, 
                             threshold=0.5, saveFile="Auto", bnet=NULL){
    
    print("Compbining results from multiple jobs and computing strenght...")
    indFileNameRes <- indFileName(dir=dir, moduleNum=moduleNum, perJob=repsPerRun, 
                                  ind=inds, typePhrase=typePhrase)

    fn <- indFileNameRes$names
    repsPerRun <- indFileNameRes$perJob
    print(paste("Looking at files with", repsPerRun, "repetition..."))
    if (!file.exists(bnInputFile)) {
        stop(paste("No file at:", bnInputFile))
    }
    bnInp <- get(load(bnInputFile)) ## bnInp
    Data <- bnInp$Data
    fn <- combinedPath(dir, fn)
    fa <- file.access(fn, 0)
    if (any(fa==-1)) {
        print("cannot access the following files, do not exist")
        print(fn[fa==-1])
        fn <- fn[fa > -1]
    }
    fa <- file.access(fn, 4)
    if (any(fa==-1)) {
        print("cannot access the following files, permission denied:")
        print(fn[fa==-1])
        fn <- fn[fa > -1]
    }
    if(length(fn)==0)
        stop("No file could be read!")
    bnCalculated <- get(load(fn[1])) ## bnCalculated
    bnet <- bnCalculated$bnet
    f1 <- bnet
    num <- repsPerRun
    absNum <- f1[, 3] * num * f1[, 4]
    for (f in fn[2:length(fn)]) {
        if (length(fn) < 2) 
            break
        if (!file.exists(f)) 
            stop(paste("No file at:", f))
        bnCalculated <- get(load(f)) ## bnCalculated
        bnet <- bnCalculated$bnet
        f1[, 3] <- bnet[, 3] + f1[, 3]
        absNum <- absNum + bnet[, 3] * repsPerRun * bnet[, 4]
        num <- num + repsPerRun
    }
    mirror <- function(x, i) {
        which(x[, 2]==x[i, 1] & x[, 1]==x[i, 2])
    }
    mirInds <- c()
    for (ind in 1:nrow(f1)) {
        mirInds[ind] <- mirror(f1, ind)
    }
    absNumMir <- absNum[mirInds]
    f1[, 3] <- f1[, 3]/length(fn)
    f1[, 4] <- absNum/(absNum + absNumMir)
    c1 <- averaged.network(strength=f1, threshold=threshold)
    c12 <- pdag2dag(c1, colnames(Data))
    sc <- bnlearn::score(c12, Data, "bde")
    scBIC <- bnlearn::score(c12, Data, "bic")
    cat(paste("**** The combined network of ", num, " models has a BDe score of ", 
              sc, " and a BIC score of ", scBIC, "\n"))
    combinedRes <- list()
    combinedRes[["bnet"]] <- f1
    combinedRes[["BN"]] <- c12
    combinedRes[["PDAG"]] <- c1
    combinedRes[["BDe"]] <- sc
    combinedRes[["BIC"]] <- scBIC
    combinedRes[["Data"]] <- Data
    combinedRes[["num"]] <- num
    combinedRes[["ind"]] <- inds
    if (doSave) {
        if (saveFile=="Auto") 
            saveFile <- get.combined.file(ind=ind, 
                                          moduleNum=moduleNum, perJob=repsPerRun, 
                                          resultPath=resdir)$combinedFile
        print(paste("Saving combined in:", saveFile))

        save.if(combinedRes, file=saveFile)
        combinedRes[["saveFile"]] <- saveFile
    }
    return(combinedRes)
}
