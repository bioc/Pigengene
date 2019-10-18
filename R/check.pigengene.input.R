check.pigengene.input <- function(Data, Labels, na.rm=FALSE, naTolerance=0.05){
    ## Performs QC on the input of one.step.pigengene() function.
    ## na.rm: If TRUE, you need to remember to update the input with the cleaned data
    ##^after running this function.
    result <- list()
    if(nrow(Data)==0)
        stop("Data has no rows!")
    if(ncol(Data)==0)
        stop("Data has no columns!")
    if(length(Labels)!=nrow(Data))
        stop("The number of rows of Data is not equal to the length of Labels!")
    if(is.null(names(Labels)))
        stop("Labels must be a named vector!")
    ## In case Lebels is a factor:
    Labels <- setNames(as.character(Labels),names(Labels))
    ## Rownames:
    if(is.null(rownames(Data)))
        stop("All input matrices must have rownames!")
    if("" %in% rownames(Data))
        stop("Rownames of the input matrices cannot be empty string ('')!")
    ## Colnames:
    if(is.null(colnames(Data)))
        stop("All input matrices must have colnames!")
    if(length(intersect(c(NA, "") , colnames(Data))) > 0)
        stop("Colnames of the input matrices cannot be NA or empty string ('')!")
    ## NAs
    Data <- check.nas(Data=Data, na.rm=na.rm, naTolerance=naTolerance)$cleaned
    if(length(intersect(rownames(Data), names(Labels)))==0)
        stop("Row of Data have no intersection with names of Labels!")
    Data <- Data[names(Labels), , drop=FALSE]
    result[["Data"]] <- Data
    result[["Labels"]] <- Labels
    return(result)
}
