check.nas <- function(
    Data, naTolerance=0.05, na.rm=TRUE)
{
    ## Checks for NAs in the Data and replaces them
    ##^with the average of the column if their frequency in the column
    ##^is not more than naTolerance.
    ## naTolerance: In the range [0, 1]. If a gene (column of Data) has more than this,
    ##^it will be removed. 
    ## na.rm: Set to FALSE to stop if there is any NA.
    result <- list()
    tooNaGenes <- c()
    nas <- is.na(Data)
    if(any(nas) & !na.rm)
        stop("Input Data has NAs!")
    tooNaCols <- which(colSums(nas)/nrow(nas) > naTolerance)
    if(any(tooNaCols)){
        Data <- Data[, -tooNaCols, drop=FALSE]
        tooNaGenes <- c(tooNaGenes, colnames(Data)[tooNaCols])
        ##genesMatched <- genesMatched[-tooNaCols]
    }
    nas <- is.na(Data)
    if(any(nas)){
        for(c1 in which(colSums(nas)>0)) ## all columns with NA
            Data[, c1][nas[, c1]] <- mean(Data[, c1], na.rm=TRUE)
    }
    replacedNaNum <- sum(nas)
    ## Warnings:
    if(length(tooNaGenes)>0)
        warning(paste(length(tooNaGenes), "genes have more than", naTolerance, 
                      "NAs and were ignored."))
    warn <- "NAs were replaced by the mean expression of the corresponding gene."
    if(replacedNaNum)
        warning(paste(replacedNaNum, warn))
    ## Output:
    result[["replacedNaNum"]] <- replacedNaNum
    result[["tooNaGenes"]] <- tooNaGenes
    result[["cleaned"]] <- Data
    result[["tooNaCols"]] <- tooNaCols
    return(result)
}##End check.nas <- function.
