save.if <- function(x1, file, verbose=1){
    result <- list()
    result[["file"]] <- file
    size <- c()
    size["inMem"] <- gdata::humanReadable(object.size(x1))
    result[["size"]] <- size
    if(is.null(file)){
        message.if(paste("Size in memory:", size["inMem"]), verbose=verbose-2)
        return(result)
    }

    ## Prepare for saving:
    nm <- deparse(substitute(x1))
    cmdToShow <- paste("save(", nm , ", file='", file, "')", sep="")
    ## R has issues with '\' but Windows needs it!
    cmd <- paste("save(", nm , ", file=file)", sep="")
    message.if(me=cmdToShow, verbose=verbose)
    assign(nm, x1)
    eval(parse(text=cmd))
    result[["size"]]["onFile"] <-  gdata::humanReadable(file.size(file))
    message.if(paste("Size in memory:", result[["size"]]["inMem"]), verbose=verbose-2)
    message.if(paste("Size on disk:", result[["size"]]["onFile"]), verbose=verbose-1)
    rm(x1)
    gc()
    return(result)
}

