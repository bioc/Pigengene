repeat.data <- function (Data, times){
    result <- list()
    if(times < 1){
        stop("times cannot be <1 !")
    }
    if (is.null(rownames(Data)))
        rownames(Data) <- 1:nrow(Data)
    repeated <- Data
    if(times >1){
        for (ind in 1:(times - 1)) {
            Datai <- Data
            rownames(Datai) <- paste(rownames(Data), ind, sep="-")
            repeated <- rbind(repeated, Datai)
        }
    }
    result[["Data"]] <- Data
    result[["repeated"]] <- repeated
    return(result)
}
