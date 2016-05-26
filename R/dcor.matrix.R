dcor.matrix <- function(Data){
    num <- ncol(Data)
    if (num > 100 * 1000) 
        stop("Not efficient for loops!")
    samples <- colnames(Data)
    if (is.null(samples)) 
        samples <- 1:num
    dcorMatrix <- matrix(NA, nrow=num, ncol=num,
                         dimnames=list(samples,samples))
    for (ind1 in 1:num)
        for (ind2 in 1:num)
            if (ind2 >=ind1) {
                if(require(energy)){
                    ed <- energy::dcor(x=Data[, ind1], y=Data[,ind2])
                } else {
                    stop("energy package is required by dcor.matrix()!")
                }
                dcorMatrix[ind1, ind2] <- ed
            }
            else {
                dcorMatrix[ind1, ind2] <- dcorMatrix[ind2, ind1]
            }
    return(dcorMatrix)
}
