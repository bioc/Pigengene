which.cluster <- function () {
    ## Determines which cluster the cod is running on.
    result <- list()
    result[["cluster"]] <- "local"
    result[["queue"]] <- NA
    result[["timeJob"]] <- Inf
    result[["numClust"]] <- 2
    result[["maxJobs"]] <- NA
    result[["onCluster"]] <- FALSE
    if (length(grep(x=system("echo $HOSTNAME", intern=TRUE), 
                    pattern="stampede")) !=0) {
        result[["cluster"]] <- "stampede"
        result[["queue"]] <- "normal"
        result[["timeJob"]] <- "47:00:00"
        result[["numClust"]] <- 15
        result[["maxJobs"]] <- 45
        result[["onCluster"]] <- TRUE
    }
    if (length(grep(x=system("echo $HOSTNAME", intern=TRUE), 
                    pattern="maverick")) !=0) {
        result[["cluster"]] <- "maverick"
        result[["queue"]] <- "gpu"
        result[["timeJob"]] <- "11:00:00"
        result[["numClust"]] <- 19
        result[["maxJobs"]] <- 18
        result[["onCluster"]] <- TRUE
    }
    return(result)
}
