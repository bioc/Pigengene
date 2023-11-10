make.filter <- function(network, epsilon, outPath=NULL){
    ##It computes a filter, which is a matrix with the
    ##same dimention as the network.
    ##If the distance between two nodes is more than epsilon,
    ##they will be disconnected in the filter, otherwise,
    ##they are connected with an edge with weight 1.
    ##network: a matrix of similarity for the network.
    ##epsilon: the threshold described above.

    result <- list()

    ##identifying a filter
    dist1 <- mean(network)/network
    filt <- dist1
    filt[dist1 <= epsilon] <- 1
    filt[dist1 > epsilon] <- 0
    diag(filt) <- 1
    
    ##Plot
    if(!is.null(outPath)){
        png(file.path(outPath, "distance_degrees.png"))
        plot(sort(log10(rowSums(dist1)/ncol(dist1))))
        dev.off()
        png(file.path(outPath, "filter_degrees.png"))
        plot(sort(log10(rowSums(filt))), main=paste("epsilon=", epsilon))
        dev.off()
    }
    
    ## Output:
    result[["epsilon"]] <- epsilon
    result[["filt"]] <- filt
    return(result)
}


    
