what.is.in <- function(x){
    if(!is.null(x)){
        load(x, verbose=TRUE)
    } else {
        return(x)
    }
}
