save.if <- function(x1, file, verbose=1){
    result <- list()
    if(is.null(file)){
        return()
    }
    nm <- deparse(substitute(x1))
    ##as.character(bquote(x))
    ##cat(nm) ## Extra --Habil
    ## assign(nm, x1)
    ##cmd <- paste("save(", nm , ", file=file)") ## Replaced with the following. --Habil. 
    cmdToShow <- paste("save(", nm , ", file='", file, "')", sep="")
    ## R has issues with '\' but Windows needs it!
    cmd <- paste("save(", nm , ", file=file)", sep="")
    message.if(me=cmdToShow, verbose=verbose)
    assign(nm, x1)
    eval(parse(text=cmd))
    size <- c()
    size["inMem"] <- gdata::humanReadable(object.size(x1))
    size["onFile"] <-  gdata::humanReadable(file.size(file))
    message.if(paste("Size in memory:", size["inMem"]), verbose=verbose-2)
    message.if(paste("Size on disk:", size["onFile"]), verbose=verbose-1)
    result[["file"]] <- file
    result[["size"]] <- size
    rm(x1)
    gc()
    return(result)
}

