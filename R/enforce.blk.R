enforce.blk <- function(bnStrength, blacklist, doWarn=TRUE)
{
    ## Removes all blacklisted arcs from the bnet.indv.
    ## Such arcs should not appear according to bnlearn specification
    ## but sometimes they do. --Habil.
    ## bnStrength: A list of objects of class bn.strength, 
    ##^ or a matrix with "from" and "to" columns.
    
    res <- list()
    sep1 <- " This-is-an-edge " ## This must not be name of any node.
    cleaned <- bnStrength
    toVector <- paste(cleaned[, 1], cleaned[, 2], sep=sep1)
    blkToVector <- paste(blacklist[, "from"], blacklist[, "to"], sep=sep1)
    inds <- which(toVector %in% blkToVector)
    if(length(inds)>0){
        cleaned <- cleaned[-inds, ]
        if(doWarn)
            warning(paste("bnStrength has", length(inds), "arcs from blacklist."))
    }
    ## Output:
    res[["inds"]] <- inds ## The list of removed indices.
    res[["cleaned"]] <- cleaned ## The cleaned list.
    return(res)
}
