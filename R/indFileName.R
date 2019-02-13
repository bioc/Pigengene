indFileName <- function (dir=NULL, moduleNum, perJob, ind, typePhrase="bnet.indv"){
    result <- list()
    if (perJob=="Auto") {
        if (is.null(dir)) 
            stop("I need the path to directory to get the number of replicates automatically!")
        where <- file.path(dir,paste("mod.", moduleNum, ".*.", "bnet.indv", 
                                     ".", ind, ".RData", sep=""))
        for (whereI in where) {
            files <- list.files(path=dirname(whereI), pattern=basename(whereI), 
                                full.names=TRUE)
            if (length(files) > 0) 
                break
        }
        if (length(files)==0) {
            print(where)
            stop("No file here to get the number of repetitions automatically!")
        }
        bnets <- get(load(files[1])) ## bnets
        bnet.indv <- bnets$bnet.indv
        repetitionShouldBe <- length(bnet.indv)
        perJob <- repetitionShouldBe
    }
    names <- paste("mod", moduleNum, perJob, typePhrase, 
                   ind, "RData", sep=".")
    result[["names"]] <- names
    result[["perJob"]] <- perJob
    ##nothing to be saved, so no name necessary
    return(result)
}
