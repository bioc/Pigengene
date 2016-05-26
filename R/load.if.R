load.if <- function(x, verbose=TRUE, ...){
    got <- NULL
    if(!is.null(x)){
        if(verbose)
            message(paste("Loading: "), x)
        y <- load(file=x, verbose=verbose, ...)
        got <- get(y)
        assign(y, got, envir=parent.frame(1))
    }
    invisible(got)
}
