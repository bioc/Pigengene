welch.pvalue <- function (cond1Values=NULL, cond2Values=NULL) {
    ## A wrapper function for t.test() and p.adjust() functions.
    ## cond1Values: matrix/dataframe containing the (expression/eigen) values for the cond1 cases
    ## cond2Values: matrix/dataframe containing the (eigen) values for the cond2 cases
    ## returns raw pvalues according to Welch's t-test, and adjusted pvalues (both by FDR and Bonferroni).
    res <- list()
    if (sum(colnames(cond1Values) !=colnames(cond2Values)) !=0) 
        stop("The expressions have not the same order in AML and MDS values!")
    ps <- vector(mode="numeric", length=ncol(cond1Values))
    for (i1 in 1:length(ps)) {
        ps[i1] <- t.test(cond1Values[, i1], cond2Values[, i1])$p.value
        ## Generalize to multi-class using stats:oneway.tes
    }
    pvals <- cbind(ps, p.adjust(ps, method="fdr"), 
                   p.adjust(ps, method="bonfer"))
    ##pvals <- cbind(colnames(cond1Values), as.data.frame(pvals))
    colnames(pvals) <- c("pValue", "FDR", "Bonferroni")
    rownames(pvals) <- colnames(cond1Values)
    res[["pvals"]] <- pvals
    return(res)
}
