dummyMaker <- function (Disease, dummies=(c(9:0) * 0.1)){
    d <- NULL
    if (!is.null(dummies)) {
        d <- matrix(rep(Disease, length(dummies), ncol=length(dummies)))
        dim(d) <- c(length(Disease), length(dummies))
        for (o in 1:length(dummies)) {
            inds <- sample(nrow(d), round(dummies[o] * nrow(d)))
            tm <- (d[inds, o] + 1)%% length(unique(Disease)) ## 2 if two-class
            d[inds, o] <- tm
        }
        colnames(d) <- paste("d", dummies, sep="")
    }
    return(d)
}
