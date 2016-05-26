make.bn.input <- function(
    moduleNum, use.Hartemink=TRUE, 
    breakNum=3, ibreaks=20, use.Disease, use.Effect, 
    saveFile, selectedFeatures=colnames(Data), 
    Data, Labels, pvalGenes=NULL, dummies=NULL, verbose=0, doPlotCor=FALSE)
{
    message.if(paste("Making BN input for module", moduleNum), verbose=verbose)
    c1 <- check.pigengene.input(Data=Data, Labels=Labels,na.rm=TRUE)
    Data <- c1$Data
    Labels <- c1$Labels
    Disease <- Labels
    Data <- Data[, selectedFeatures]
    rownamesData <- rownames(Data) ## discretize() may change them!
    if (moduleNum=="E" & doPlotCor){
        corPlotted <- draw.cor.cond(Data=Data, Labels=Labels, verbose=verbose-1,
                                    savePath=dirname(saveFile))
    }
    if (!use.Hartemink) {
        Data <- discretize(data=as.data.frame(Data), 
                           breaks=breakNum, method="interval")
    }
    else {
        Data <- discretize(data=as.data.frame(Data), 
                           breaks=breakNum, method="hartemink", ibreaks=ibreaks)
    }
    d <- NULL
    if(!is.null(dummies)){
        if(length(unique(Disease))!=2){
            warning("Dummies for multiclass case are experimental. There are more than two disease classes.")
            ##dummyMaker now works with more than 2 classes, but lets say is 'experimental'.
        }
        d <- dummyMaker(Disease=as.numeric(as.factor(Disease)), dummies=dummies)
    }
    if (use.Disease) {
        ##d <- dummyMaker(Disease, dummies)
        d <- as.data.frame(cbind(d, Disease))
        d$Disease <- as.factor(d$Disease)##just to be sure
    }
    if (use.Effect) {
        Effect <- Disease
        d <- as.data.frame(cbind(Effect, d))
        d$Effect <- as.factor(d$Effect)##just to be sure
    }
    if (use.Disease | use.Effect) {
        Data <- cbind(Data, discretize(as.data.frame(d), method="interval", 
                                       breaks=length(unique(Disease)))) ## Most likely binary.
    }
    blk <- c()
    if (use.Disease) {
        blk <- rbind(blk, tiers2blacklist(list("Disease", setdiff(colnames(Data), 
                                                                  "Disease"))))
    }
    if (use.Effect) {
        blk <- rbind(blk, tiers2blacklist(list(setdiff(colnames(Data), 
                                                       "Effect"), "Effect")))
    }
    rownames(Data) <- rownamesData
    result=list()
    result[['Data']] <- Data
    result[['blacklist']] <- blk 
    bnInp <- result
    save.if(bnInp, file=saveFile, verbose=verbose)
    message.if(paste("BN input file was saved in:", saveFile), verbose=verbose)
    return(bnInp)
}
