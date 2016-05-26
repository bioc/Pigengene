getBN <- function (x, nodelist, verbose=0){
    ## x: c(<The file containing bnets>, <repeatition>)
    bnets <- get(load(x[1])) ## bnets
    bnet.indv <- bnets$bnet.indv
    message.if(x[1], verbose=verbose)
    rep <- as.numeric(x[2])
    emptBN <- empty.graph(nodes=nodelist)
    candy <- emptBN
    bnlearn::arcs(candy) <- bnet.indv[[rep]]
    return(candy)
    ## absolutely no point in a list
}
