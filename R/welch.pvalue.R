welch.pvalue <- function (Data, Labels){ ##cond1Values=NULL, cond2Values=NULL) {
    ## A wrapper function for oneway.test (which generalizes t.test() to multiple groups)
    ## and p.adjust() functions.
    ## Data: matrix/dataframe. Each colum contains the (expression/eigen) values for a gene.
    ## returns raw pvalues according to Welch's t-test, and adjusted pvalues (both by FDR and Bonferroni).
    res <- list()
    if (sum(rownames(Data) !=names(Labels)) !=0) 
        stop("The rows of Data should be the same as the names of Labels with the same order!")
    if(is.null(colnames(Data)))
        colnames(Data) <- paste("V",1:ncol(Data),sep="")
    ps <- vector(mode="numeric", length=ncol(Data))
    for (i1 in 1:length(ps)) {
        ##ps[i1] <- t.test(cond1Values[, i1], cond2Values[, i1])$p.value
        ## Generalize to multi-class using stats:oneway.tes
        f1 <- as.formula(paste(colnames(Data)[i1],"~ Labels",sep=" "))
        ps[i1] <- oneway.test(data=Data,formula=f1)$p.value
    }
    pvals <- cbind(ps, p.adjust(ps, method="fdr"), 
                   p.adjust(ps, method="bonfer"))
    ##pvals <- cbind(colnames(cond1Values), as.data.frame(pvals))
    colnames(pvals) <- c("pValue", "FDR", "Bonferroni")
    rownames(pvals) <- colnames(Data)
    res[["pvals"]] <- pvals
    return(res)
}
