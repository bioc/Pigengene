combinedPath <- function (dir=NULL, fn)
{
    if(is.null(fn)){
        return(fn)
    }
    if (is.null(dir)) {
        dir <- getwd()
    }
    dir <- normalizePath(dir)
    fnc <- file.path(dir, fn)
    return(fnc)
    ## There is absolutely no point in returning a list.
}
