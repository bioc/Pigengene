get.genes <- function(
    c5Tree=NULL, pigengene=NULL, queue=NULL, modules=NULL, pos=0, enhance=TRUE)
{
    ## needs either: [net and queue] OR [c5Tree and pigengene]
    ## returns all genes that are left after shrinking for 'pos' steps
    ## if 'enhance', makes sure that the output contains at least two genes from each used module .
    
    if(!is.null(c5Tree) & ! is.null(pigengene) ){
        feats <- get.used.features(c5Tree)
        modules <- pigengene$orderedModules[colnames(pigengene$Data)]
        queue <- make.membership.queue(feats, pigengene, modules=modules) 
    }
    if(is.null(modules) | is.null(queue) ){
        stop("Provide (modules and queue) or (c5Tree and  pigengene)!")
    }
    
    if(enhance){
        an <- modules[names(queue)]
        ## lets take a note of the two genes with the highest membership in each module.
        last2s <- c()
        for(h1 in unique(an)){
            last2s <- c(last2s, tail(names(queue)[which(modules[names(queue)]==h1)], 2))
        } 
        a1p <- names(queue)[-(1:pos)]
        amiss <- setdiff(last2s, a1p)
        mylis <- c(amiss, a1p)
    }else{
        mylis <- names(queue)[-(1:pos)]}
    return(mylis)
}


