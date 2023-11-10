module.heatmap <- function(c5Tree=NULL, pigengene, mes=NULL, saveDir, testD=NULL, testL=NULL, 
                           pos=0, verbose=0, doAddEigengene=TRUE,  scalePngs=1, ...){
    ## Takes a decision tree and a pigengene as input
    ## and creates one heatmap for each module
    ## pos > 0 removes genes similar to compact.tree behavior
    result <- list()
    result[["call"]] <- match.call()
    result[["saveDir"]] <- saveDir
    trainDir <- combinedPath(saveDir, 'train')
    testDir <- combinedPath(saveDir, 'test')
    allDir <- combinedPath(saveDir, 'all')

    ## Expression data:
    D1 <- pigengene$Data
    
    ## QC:
    if(!is.null(c5Tree)){
        message.if("Plotting the heatmaps for the nodes of the tree, mes ignored.",
                   verbose=verbose)
        feats <- get.used.features(c5Tree=c5Tree)
        plotDir <- trainDir
    } else { ## no c5Tree,
        feats <- colnames(pigengene$eigengenes)
        if(is.null(mes)){
            message.if("Plotting the heatmaps for all modules...", verbose=verbose)
        } else {
            feats <- intersect(feats, mes)
        }
        plotDir <- allDir
    }
    message.if("Features:", verbose=verbose-1)
    message.if(paste(feats, collapse=", "), verbose=verbose-1)
    ## QC:
    if(length(feats)==0)
        stop("No features to plot!")
    modules0 <- pigengene$orderedModules
    ## It is possible that a gene is not present in the D1 (pigengene$Data), but be in the modules,
    ##^ for example, when two datasets are used and the gene is NA in one of the datasets.
    ## We cannot include such genes in the heatmap.
    modules <- modules0[intersect(names(modules0), colnames(D1))]
    if(length(modules0) > length(modules)){
        warning(paste("Data are not available for", length(modules0)-length(modules), "genes in the modules."))
    }
    ## If we are only interested in the expression of the genes in the compact Tree:
    if(pos>0 & !is.null(c5Tree)){
        whatMatters <- get.genes(pos=pos, c5Tree=c5Tree, enhance=TRUE, pigengene=pigengene)
        modules[-match(whatMatters, names(modules))] <- -1
    }
    ## drop constant and NA columns --Amir
    rmnas <- function(dat){
        sds <- sapply(as.data.frame(dat), sd)
        if(any(sds<1E-7 | is.na(sds))){
            dat <- dat[, -which(sds<1E-7 |is.na(sds)), drop=FALSE]
        }
        return(dat)
    }
    minmax <- function(x){
        if(min(x)==max(x)){ret <- rep(0, length(x))}
        if(min(x)<max(x)){
            ret <- ((x-min(x))/(max(x)-min(x)))}
        return(ret)}
    mkpngs <- function(f1, Data, saveDir1, anR){
        if(ncol(Data) - doAddEigengene <2){
            errorMessage <- paste("I need at least 2 columns (genes) to plot a heatmap.",
                                  f1, "skipped.")
            message.if(errorMessage, verbose=verbose-3)
            return(NULL)
        }
        dir.create(path=saveDir1, recursive=TRUE, showWarnings=FALSE)
        message.if(paste("Plotting heatmaps for",f1,"with",ncol(Data),"columns in:", saveDir1),
                   verbose=verbose)
        legwidth=16*max(nchar(as.character(anR)))
        ##legwidth=16*max(nchar(as.character(pigengene$annotation[, 1])))
        if(scalePngs!=1){
            wdt <- 220+(scalePngs*ncol(Data))+legwidth  ## scalePngs used to be 7 here.
            hgt <- 210+(scalePngs*nrow(Data)) ##scalePngs used to be 7.8 here.
            fontsizeH <- scalePngs*20
        } else {
            wdt <- 4*480
            hgt <- 2*480
            fontsizeH <- 20
        }
        fn1 <- combinedPath(dir=saveDir1, fn=paste0(f1, '.png'))
        message.if(paste("Plotting:", fn1), verbose=verbose-1)      
        png(fn1, width=wdt, height=hgt)
        pheatmap.type(Data=Data, annRow=anR, type='type', fontsize=fontsizeH, show_rownames = FALSE, ...)
        dev.off()
        fn2 <- combinedPath(dir=saveDir1, fn=paste0(f1, "_scaled.png"))
        message.if(paste("Plotting:", fn2), verbose=verbose-1)
        png(fn2, width=wdt, height=hgt)
        pheatmap.type(Data=scale(Data), annRow=anR, type='type',fontsize=fontsizeH, ...)
        dev.off()
        fn3 <- combinedPath(dir=saveDir1, fn=paste0(f1, "_scaledMinMax.png"))
        message.if(paste("Plotting:", fn3), verbose=verbose-1)
        png(fn3, width=wdt, height=hgt)
        pheatmap.type(Data=apply(scale(Data), FUN=minmax, MARGIN=2), annRow=anR, type='type',
                      fontsize=fontsizeH, ...)
        dev.off()
    }
    
    anR <- pigengene$annotation
    colnames(anR) <- 'type'
    ##Check
    if(is.null(rownames(anR)))
        stop("Not a valid pigengene object because 'annotation' does not have row names!")
    if(is.null(rownames(D1)))
        stop("Not a valid pigengene object because 'Data' does not have row names!")
    if(any(sort(rownames(anR))!=sort(rownames(D1))))
        stop("Not a valid pigengene object because row names differ in 'Data' and annotations!")
    anR <- anR[rownames(D1),,drop=FALSE]

    ## for testD
    D1test <- NULL
    anRtest <- NULL
    if(!is.null(testD) & !is.null(testL)){
        D1test <- testD
        anRtest <- as.data.frame(testL)
        colnames(anRtest) <- 'type'
        anRtest <- anRtest[rownames(D1test),,drop=FALSE]
    }
    
    ## combined expression plot of the genes in relevant modules
    if(!is.null(c5Tree)){
        ## Subsetting only the genes in the modules:
        Data1 <- rmnas(D1[, names(modules)[which(modules %in%gsub("ME", "", feats))]])
        mkpngs(f1="combined", Data=Data1, saveDir1=plotDir, anR=anR)
        if(!is.null(anRtest)){
            testdata <- rmnas(D1test[, match(colnames(Data1), colnames(D1test))])
            mkpngs(f1="combined", Data=testdata, saveDir1=testDir, anR=anRtest)
        }
    }

    ## Per module plots
    for(f1 in feats){
        message.if(paste("Module:", f1), verbose=verbose-1)
        Data1 <- rmnas(D1[, names(modules)[which(modules==gsub("ME", "", f1))], drop=FALSE])
        if(doAddEigengene)
            Data1 <- cbind(Data1, pigengene$eigengenes[rownames(Data1),f1,drop=FALSE])
        mkpngs(f1=f1, Data=Data1, saveDir1=plotDir, anR=anR)
        if(!is.null(anRtest)){
            testdata <- rmnas(D1test[, match(colnames(Data1), colnames(D1test))])
            mkpngs(f1=f1, Data=testdata, saveDir1=testDir, anR=anRtest)
        }
    }
    invisible(result)
}
