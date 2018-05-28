save.if <- function(x1, file, verbose=1){
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
    message.if(me=paste(cmdToShow, '\n'), verbose=verbose)
    assign(nm, x1)
    eval(parse(text=cmd))
}

