get.used.features <- function(c5Tree){
    if(!inherits(c5Tree,"C5.0"))
        stop("The class of c5Tree argument must be 'C5.0' !")
    
    v1 <- unlist(strsplit(c5Tree$tree, 'att'))[-1]
    v2 <- c()
    for(i in 1:length(v1)){
        v2 <- c(v2, unlist(strsplit(v1[i], " "))[1])
    }
    
    v3 <- gsub("=", "", v2)
    v4 <- gsub('"', "", v3)
    ##need unique here, in case the tree reuses a feature 
    return(unique(v4))
}
