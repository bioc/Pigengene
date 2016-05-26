draw.bnS <- function(
    inputType="ENTREZIDat", edgeColor="blue", DiseaseCol="darkgreen", 
    DiseaseFill="green", doShowBDe=FALSE, verbose=0, 
    consensusFile, consensusThresh, moduleNum="E", moduleNamesFile=NULL) 
{
    ## moduleNamesFile: a csv file that maps modules names ("ME28")
    ## to more meaningful names (e.g. "HOX"). --Habil.
    result <- list()
    plots <- list()
    mainText <- NULL
    bde <- c()
    result[["moduleNamesFile"]] <- moduleNamesFile
    message.if("Plotting graphs...", verbose=verbose)
    if (moduleNum=="E" & !is.null(moduleNamesFile)) {
        if (!file.exists(moduleNamesFile)){
            warning(paste("moduleNamesFile does not exist at:", moduleNamesFile, "Ignored."))
            moduleNamesFile <- NULL
        }
    }
    else {
        moduleNamesFile <- NULL
    }
    files <- get.consensus.files(consensusFile=consensusFile, 
                                 threshVector=consensusThresh)##StrengthsThres
    for (fileI in files) {
        consensusRes <- get(load(fileI, verbose=verbose>0)) ## consensusRes
        bde[fileI] <- consensusRes$BDe
        if (doShowBDe) 
            mainText <- paste("BDe=", round(bde[fileI]))
        plotFile <- gsub(fileI, pattern="\\.RData$", replacement=".png")
        plots[[fileI]] <- draw.bn(BN=consensusRes$BN, verbose=verbose-1, 
                                  plotFile=plotFile, 
                                  moduleNamesFile=moduleNamesFile, mainText=mainText)
    }
    result[["files"]] <- files
    result[["plots"]] <- plots
    result[["bde"]] <- bde
    return(result)
}
