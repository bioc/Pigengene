make.decision.tree <- function(
    pigengene, Data=pigengene$Data,
    Labels=structure(pigengene$annotation[rownames(pigengene$eigengenes),1],
        names=rownames(pigengene$eigengenes)),
    testD=NULL, testL=NULL, 
    selectedFeatures=NULL, saveDir='C5Trees', minPerLeaf=NULL, 
    useMod0=FALSE, costRatio=1, toCompact=NULL, noise=0, noiseRepNum=10, 
    doHeat=TRUE, verbose=0)
{
    ##^Make a decision tree.
    ##pigengene: from compute.pigengene, the $eigengenes part are the potential
    ##^predictors. The $orderedModules part is needed by compact.tree
    ## selectedFeatures : if NULL, all partial Eigengenes are treated as potential
    ##^predictors. Otherwise, only the ones listed in selectedFeatures are considered for
    ##^building the tree(e.g. the children of Disease in the consensus BN).
    ##Data(1, 2), testD(1, 2) the actual gene expression profiles, 
    ##^used for projection (compacting the tree), and used for creating the response
    ##^vectors
    ## minPerLeaf : The desired minimal number of samples in each leaf node.
    ##^By default (set to NULL) we try out all values between 2% and 10% of the total
    ##^number of samples.
    ## toCompact: A single value in minPerLeaf that is used for compacting the tree.
    ##^ Set to FALSE not to compact.
    ## noise: A value in [0, 1]. If not 0, upto this portion of the test expression
    ## will be replaced by Gaussian noise to estimate sensitivity to noise.
    ## noiseRepNum: Number of repeats for computing the average accuracy over noisy data.
    ## took out: Data1=NULL, Data2=NULL , testD1=NULL,  testD2=NULL, 
    message.if("Making decision trees...", verbose=verbose)
    Data <- Data[names(Labels), ]
    ##Data <- Data[rownames(pigengene$eigengenes), ]
    ## QC:
    if(class(pigengene)!="pigengene")
        stop("pigengene must be an object of class pigengene!")
    c1 <- check.pigengene.input(Data=Data, Labels=Labels, na.rm=TRUE)
    Data <- c1$Data
    Labels <- c1$Labels
    ##
    c5Res <- list()
    c5Res[["call"]] <- match.call()
    dir.create(path=saveDir, recursive=TRUE, showWarnings=FALSE)
    if(is.null(testD)){
        ## test is same as train
        testD <- Data
        testL <- Labels 
    }
    if(is.null(minPerLeaf)){
        ltt <- length(Labels)
        ## For really unbalanced datasets, we need to reduce the maximum of minNodesperLeafs
        ltmin <- round(0.7*min(table(Labels)))
        minPerLeaf=unique(ceiling(0.02*ltt):(min(ltmin, ceiling(0.1*ltt))))
    }
    m1 <- paste(minPerLeaf, collapse=", ")
    message.if(me=paste("minPerLeaf:", m1), verbose=verbose)
    if(is.null(selectedFeatures)){
        selectedFeatures <- colnames(pigengene$eigengenes)
    }
    if(!useMod0)
        selectedFeatures <- setdiff(selectedFeatures, "ME0")
    if(length(unique(Labels))>2 & costRatio!=1){
        warning("costRatio for multiClass case not supported yet, reverting to costRatio=1.")
    }

    inpD <- as.data.frame(cbind(pigengene$eigengenes[names(Labels), selectedFeatures], Labels))
    assign('inpD', inpD , envir=parent.env(parent.frame()), inherits=TRUE)
    colnames(inpD) <- c(selectedFeatures, 'Labels')
    ## Costs:
    mycosts <- matrix(1, nrow=length(unique(Labels)), ncol=length(unique(Labels)))
    diag(mycosts)=0
    if(length(unique(Labels))==2 & costRatio!=1){  
        mycosts <- matrix(c(0, costRatio, 1, 0), nrow=2, ncol=2)
    }
    colnames(mycosts) <- levels(as.factor(inpD$Labels))
    rownames(mycosts) <- levels(as.factor(inpD$Labels))
    if(verbose>0){
        message("costs:")
        print(mycosts)
    }
    c5TreeS <- list()
    fitD <- c()
    ##for plotting: 
    pngplot <- function(c5Tr, h1, saveDir=NULL){
        lbls <- as.character(inpD$Labels)
        mycols <- rgb(expand.grid(c(0:1), (0:1), (0:1)))[c(2, 5, 3, 1, 4, 6:7)]
        errQs <- table(predict(c5Tr, inpD), lbls)
        errs <- colSums(errQs)-diag(errQs)
        Acc <- 1-round((sum(errs)/sum(errQs))*100)/100## accuracy
        adf <- (log2(log2(nchar(c5Tr$output))))
        extr <- adf*sum(nchar(unique(lbls)))*max(1, log10(1+length(grep(" freq", unlist(strsplit(c5Tr$tree, "=")))))^3)
        pngfn <- combinedPath(dir=saveDir, fn=paste('C5tree', h1, '.png', sep=''))
        myinfo <- paste("minCases:", h1 , "Acc:", Acc , "Misclassified:" , paste(errs, names(errs), sep=' ', collapse='...'))
        png(pngfn, height=max(1, adf*0.6)*480, width=extr+(max(1, adf*0.8)*480))
        plot.new()
        temp <- legend("topleft", legend="", inset=-0.005, bty='n')
        plot(c5Tr, gp=grid::gpar(col ="blue"), legend.=temp, 
             main=myinfo, tp_args=list("fill" = mycols))
        temp <- legend("topleft", legend = unique(lbls)[order(unique(lbls))], 
                       text.col=mycols, inset=-0.005, bty='n', ncol=ceiling(length(unique(lbls))/4))
        dev.off()
        return(myinfo)
    }
    
    ##We build several trees at different thresholds for minimal number of cases per leaf, 
    ##^ either using all pigenegenes as input, or only a selected few.
    for(h1 in minPerLeaf){
        c5cnt <- C5.0Control(minCases=h1)
        ## Build a tree that models the Labels as a function of the pigenegenes, 
        ##^ there should be at least h1 cases in each leaf, do not boost.
        c5Tr <- C5.0(as.factor(Labels)~., data=inpD, control=c5cnt, trials=1, costs=mycosts)
        ##fn <- combinedPath(dir=saveDir, fn=paste('C5tree', h1, '.RData', sep=''))
        ##save(inpD, c5Tr, file=fn)
        browser()
        fitD <- cbind(fitD, (unlist(get.fitted.leaf(c5Tr, inpD))))
        colnames(fitD)[ncol(fitD)] <- h1  
        c5Tr$info <- pngplot(c5Tr, h1, saveDir)
        c5TreeS[[as.character(h1)]] <- c5Tr ## Habil named with characters.
    }
    ## Compacting:
    if(is.null(toCompact)){
        ## ignore trees that do not have any leafs, find most general proper tree
        atmp <- max(minPerLeaf)
        while(is.na(get.used.features(c5TreeS[[as.character(atmp)]])) & atmp >min(minPerLeaf)){
            atmp=atmp-1
        }
        if(!is.na(get.used.features(c5TreeS[[as.character(atmp)]]))){
            toCompact <- atmp}
        else{ ##no proper trees found
            toCompact <- FALSE }
    }
    message.if(paste("toCompact:", toCompact), verbose=verbose)
    if(toCompact){
        tc5 <- c5TreeS[[as.character(toCompact)]]
        compactC5Tr <- compact.tree(c5Tree=tc5, pigengene=pigengene, 
                                    Data=Data, Labels=Labels, 
                                    testD=testD, testL=testL, 
                                    saveDir=saveDir, verbose=verbose-1)
        c5Res[['compacted']] <- compactC5Tr
        if(doHeat){
            heatCompact <- module.heatmap(
                c5Tree=tc5, pigengene=pigengene, 
                saveDir=combinedPath(saveDir, "heatmaps_compact"), 
                testD=testD, testLabels=testL, 
                pos=compactC5Tr$pos, verbose=verbose-1)
            c5Res[['heatCompact']] <- heatCompact
        }
    } else ## Do not compact.
        tc5 <- c5TreeS[[as.character(max(minPerLeaf))]]
    ## ^^ why is this here?! --Amir 
    ## Noise:
    if(noise!=0){
        noiseDir <- combinedPath(saveDir, "noise")
        noisy <- noise.analysis(c5tree=tc5, pigengene=pigengene, 
                                Data=Data, verbose=verbose, 
                                testD=testD, testL=testL, 
                                saveDir=noiseDir, noise=noise, repNum=noiseRepNum)
        c5Res[['noisy']] <- noisy
    }
    ## Plot heatmaps:
    if(doHeat){
        heat <- module.heatmap(
            c5Tree=tc5, pigengene=pigengene, 
            saveDir=combinedPath(saveDir, "heatmaps"), 
            testD=testD, testLabels=testL, pos=0, 
            verbose=verbose-1)
        c5Res[['heat']] <- heat
    }
    ## Output:
    c5Res[['minPerLeaf']] <- minPerLeaf
    c5Res[['c5Trees']] <- c5TreeS
    c5Res[['leafLocs']] <- fitD
    c5Res[['toCompact']] <- toCompact
    c5Res[['costs']] <- mycosts
    c5Res[['saveDir']] <- saveDir
    return(invisible(c5Res))
}
