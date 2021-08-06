type2char <- function(type1){
    ## Identify the class of input.
    if(inherits(type1, "character")){
        if(type1=="Human"){
            char1 <- "org.Hs.eg.db"
        } else if(type1=="Mouse"){
             char1 <- "org.Mm.eg.db"
        } else{
            stop("A character datatype other than 'Human' or 'Mouse' is not supported!")
        }
    } else {
        char1 <- type1$packageName
    }
    return(char1)
}
