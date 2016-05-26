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
    ##^H.
    if(FALSE){ 
        str <- dir
        if (substr(str, (nchar(str) + 1) - 1, nchar(str))=="/") {
            fnc <- paste(str, fn, sep="")
        }
        else {
            fnc <- paste(str, fn, sep="/")
        }
    }
    return(fnc)
    ## there is absolutely no point in returning a list.
}
