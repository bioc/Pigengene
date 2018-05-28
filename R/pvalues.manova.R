pvalues.manova <- function(Data, Labels){
    res <- list()
    res[["call"]] <- match.call()
    Data <- as.matrix(Data)
    fit <- manova(Data~as.factor(Labels))
    smry <- unlist(summary.aov(fit))
    ps <- smry[grep("F)1", names(smry))]
    pvals <- cbind(ps, p.adjust(ps, method="fdr"), 
                   p.adjust(ps, method="bonfer"))
    colnames(pvals) <- c("pValue", "FDR", "Bonferroni")
    rownames(pvals) <- colnames(Data)
    res[["pvals"]] <- pvals
    res[["manovaFit"]] <- fit
    return(res)
}
