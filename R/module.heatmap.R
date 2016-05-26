module.heatmap <- function(c5Tree, pigengene, saveDir, testD=NULL, testL=NULL, 
                           pos=0, verbose=0, ... ){
    ## Takes a decision tree and a pigengene as input
    ## and creates one heatmap for each module
    ## pos > 0 removes genes similar to compact.tree behaviour
    result <- list()
    result[["call"]] <- match.call()
    result[["saveDir"]] <- saveDir
    trainDir <- combinedPath(saveDir, 'train')
    testDir <- combinedPath(saveDir, 'test')
    if(!is.null(c5Tree)){
        feats <- get.used.features(c5Tree=c5Tree)
    } else {
        feats <- colnames(pigengene$eigengenes)
    }
    modules <- pigengene$orderedModules
    ## If we are only interested in the expression of the genes in the compact Tree:
    if(pos>0 & !is.null(c5Tree)){
        whatMatters <- get.genes(pos=pos, c5Tree=c5Tree, enhance=TRUE, pigengene=pigengene)
        modules[-match(whatMatters, names(modules))] <- -1
    }
    ## drop constant and NA columns --Amir
    rmnas <- function(dat){
        sds <- sapply(as.data.frame(dat), sd)
        if(any(sds<1E-7 | is.na(sds))){
            dat <- dat[, -which(sds<1E-7 |is.na(sds))]
        }
        return(dat)
    }
    minmax <- function(x){
        if(min(x)==max(x)){ret <- rep(0, length(x))}
        if(min(x)<max(x)){
            ret <- ((x-min(x))/(max(x)-min(x)))}
        return(ret)}
    mkpngs <- function(f1, data, saveDir1, anR){
        dir.create(path=saveDir1, recursive=TRUE, showWarnings=FALSE)
        message.if(paste("Ploting heatmaps in:", saveDir1), verbose=verbose)
        fn <- paste(f1, '.png', sep='') ## Habil.
        ##fn=paste('heat', f1, '.png', sep='')
        legwidth=16*max(nchar(as.character(pigengene$annotation[, 1])))
        wdt <- 220+(7*ncol(data))+legwidth
        hgt <- 210+(7.8*nrow(data))
        png(combinedPath(dir=saveDir1, fn=fn), width=wdt, height=hgt)
        pheatmap.type(Data=data, annRow=anR, type='type', fontsize_col=7, show_rownames = FALSE, ...)
        dev.off()
        fn2 <- paste("scaled_", fn, sep='')
        png(combinedPath(dir=saveDir1, fn=fn2), width=wdt, height=hgt)
        pheatmap.type(Data=scale(data), annRow=anR, type='type', ...)
        dev.off()
        fn3 <- paste("scaledMinMax_", fn, sep='')
        png(combinedPath(dir=saveDir1, fn=fn3), width=wdt, height=hgt)
        pheatmap.type(Data=apply(scale(data), FUN=minmax, MARGIN=2), annRow=anR, type='type', ...)
        dev.off()
    }
    
    D1 <- pigengene$Data
    anR <- pigengene$annotation
    colnames(anR) <- 'type'
    rownames(anR) <- rownames(D1)
    ## for testD
    D1test <- NULL
    anRtest <- NULL
    if(!is.null(testD)&!is.null(testL)){
        D1test <- testD
        anRtest <- as.data.frame(testL)
        colnames(anRtest) <- 'type'
        rownames(anRtest) <- rownames(D1test) 
    }
    ## combined expression plot of the genes in relevant modules 
    if(!is.null(c5Tree)){
        data <- rmnas(D1[, names(modules)[which(modules %in%gsub("ME", "", feats))]]) 
        mkpngs(f1="combined", data=data, saveDir1=trainDir, anR=anR)
        if(!is.null(anRtest)){
            testdata <- rmnas(D1test[, match(colnames(data), colnames(D1test))])
            mkpngs(f1="combined", data=testdata, saveDir1=testDir, anR=anRtest)
        }
    }
    ## Per module plots
    for(f1 in feats){
        ##long lines :)
        data <- rmnas(D1[, names(modules)[which(modules==gsub("ME", "", f1))]])
        mkpngs(f1=f1, data=data, saveDir1=trainDir, anR=anR)
        if(!is.null(anRtest)){
            testdata <- rmnas(D1test[, match(colnames(data), colnames(D1test))])
            mkpngs(f1=f1, data=testdata, saveDir1=testDir, anR=anRtest)
        }
    }
    invisible(result)
}
