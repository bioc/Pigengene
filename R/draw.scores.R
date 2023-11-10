draw.scores <- function(candidates, maxUpto=1500, candlist=NULL, verbose=0, savePath){
    result <- list()
    scores <- as.numeric(candidates[, "Score"])
    num <- length(scores)
    png(file.path(savePath, paste("scores.hist", num, ".png", sep="")))
    hist(scores, main=paste("Histogram of", num, "scores\n Max=", 
                     max(scores)))
    dev.off()
    png(file.path(savePath, paste("scores.improvement", num, ".png", sep="")))
    plotted <- draw.improvement(scores=scores[1:min(length(scores), 
                                    maxUpto)])
    dev.off()
    message.if(paste("Two score plots were saved in:", savePath), verbose=verbose)
    result[["maxes"]] <- plotted$maxes
    result[["scores"]] <- scores
    return(result)
}
