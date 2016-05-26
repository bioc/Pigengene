sample.network <- function(
    nodes, blacklist=NULL, revBlkRate=0, doShuffle=TRUE, 
    prob=min(2 * 2/(length(nodes)-1), 1), verbose=0)
{
    ## Returns a BN, to be used as a starting point in hc.-- Amir
    ## blk is a blacklist (matrix/dataframe with from, to columns)
    ## revBlkRate is the fraction (between 0 and 1) of reversed black listed arcs 
    ##^ that should be added in reverse direction. 
    ## We are going to uase the 'ordered' method, so lets first permutate the nodes
    ##^ randomly to make sure that the order is arbitrary
    ## prob: The probability of an edge being in the graph.
    ##^The bnlearn default is 2/(length(nodes)-1), 
    ## which gives the expected number of edges to be equal to the # of nodes.
    message.if("sampling.network...", verbose=verbose)
    if(is.null(nodes))
        stop("nodes cannot be NULL!")
    if(prob<0 | 1 <prob)
        stop("prob must be a numeric value in [0, 1] !")
    message.if(paste("doShuffle:", doShuffle), verbose=verbose-1)
    if (doShuffle){
        message.if("Nodes ordering before shuffling:", verbose=verbose-2)
        message.if(nodes, verbose=verbose-2)
        nodes <- sample(nodes, length(nodes), replace=FALSE)
        if("Disease" %in% nodes) ## Move Disease to the begining.
            nodes <- c("Disease", nodes[-which(nodes=="Disease")])
    }
    F1 <- random.graph(nodes, method="ordered", num=1, prob=prob)
    message.if("Nodes ordering matters for random graph:", verbose=verbose-2)
    message.if(nodes, verbose=verbose-2)
    if (revBlkRate > 0 & !is.null(blacklist)) {
        n <- round(revBlkRate * nrow(blacklist))
        inds <- sample(nrow(blacklist), n)
        revs <- blacklist[inds, c("to", "from")]
        colnames(rev) <- c("from", "to")
        f1rev <- rbind(arcs(F1), rev)
        arcs(F1) <- f1rev
    }
    if (!is.null(blacklist)) {
        f1b <- rbind(arcs(F1), blacklist)
        if (any(duplicated(f1b))) {
            toremove <- which(base::duplicated(f1b, fromLast=TRUE))
            arcs(F1) <- arcs(F1)[-toremove, ]
        }
    }
    return(F1)
}
