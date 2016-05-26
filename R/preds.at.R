preds.at <- function(c5Tree, pigengene, pos=0, Data){
    ## returns vector of prediction of c5Tree on Data after 'pos' steps of shrinking
    pi2 <- pigengene
    feats <- get.used.features(c5Tree)
    modules <- pigengene$orderedModules[colnames(pigengene$Data)]
    queue <- make.membership.queue(feats, pigengene, modules=modules)

    genes <- get.genes(queue=queue, pos=pos, enhance=TRUE, modules=modules)
    pi2$orderedModules[setdiff(names(queue), genes)] <- -1
    ##pi2$orderedModules[setdiff(names(queue), get.genes(queue=queue, pos=pos, enhance=TRUE, modules, c5Tree))] <- -1
    ##queue, pos, enhance=TRUE, modules=NULL, c5Tree=NULL
    rr2 <- project.eigen(Data=Data, pigengene=pi2, verbose=1, ignoreModules= -1 )
    res <- list()
    res[["predictions"]] <- predict(c5Tree, as.data.frame(rr2$projected))
    res[["eigengenes"]] <- as.data.frame(rr2$projected)
    return(res)
}
