draw.improvement <- function(scores, ...) {
    result <- list()
    maxes <- c()
    for (ind in 1:length(scores)) {
        maxes[ind] <- max(scores[1:ind])
    }
    plot(maxes, main=paste("Improvement of scores"), type="l", ...)
    result[["maxes"]] <- maxes
    return(result)
}
