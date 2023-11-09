apply.filter <- function(gamma, filt, Data, doNormalize=FALSE){
    ## It applies a filter on the Data, which is a matrix with the
    ## same dimention as the Data by using this formula:
    ## (gamma*Data*filt)+(1-gamma)*Data
    ## gamma: This value is used in above formula which was set in runall.R.
    ## filt: the filter matrix which is the output of identify.filter().
    ## Data: the gene experission matrix with samples on rows.
    ## Output: filtered, the filtered matrix.
    
    result <- list()
    if(inherits(Data, "list")){##Data is a list
        filtered <- list()
        print(paste("gamma", gamma))
        for(i1 in 1:length(Data)){
            d1 <- Data[[i1]]
            f1 <- apply.filter(gamma=gamma, filt=filt, Data=d1)
            filtered[[i1]] <- f1$filtered
            
        }
    } else {##Data is a matrix
        ## QC	
        if(any(!colnames(Data) %in% colnames(filt)))
            stop("filt and Data must have the same number of columns!")
        Data <- as.matrix(Data)
        ##Applying a filter on the Data
        meanData <- colMeans(Data)
        DataSds <- colSds(Data)
        sData <- scale(Data)
        filtN <- filt[colnames(sData), colnames(sData)]
        ##Normalize the filter by the degree in the graph
        if(doNormalize){
            filtN <- filtN / rep(rowSums(filtN), each=nrow(filtN))
        }
        filtered <- (gamma*sData %*% filtN) + ((1-gamma)*sData)
        filtered <- t(t(filtered) * DataSds)
        filtered <- t(t(filtered) + meanData)
    }
    ## Output:
    result[["filtered"]] <- filtered
    return(result)
}
