indFileName <- function(dir=NULL, moduleNum, perJob, ind, typePhrase="bnet.indv", verbose=0){
    result <- list()
    message.if("indFileName", verbose=verbose)
    if(perJob=="Auto"){
        message.if(dir, verbose=verbose-1)
        if(is.null(dir)) 
            stop("I need the path to directory to get the number of replicates automatically!")
        where <- file.path(dir, paste("mod.", moduleNum, ".*.", "bnet.indv", 
                                     ".", ind, ".RData", sep=""))
        for(whereI in where) {
            files <- list.files(path=dirname(whereI), pattern=basename(whereI), 
                                full.names=TRUE)
            if(length(files) > 0) 
                break
        }
        if(length(files)==0) {
            m1 <- "No file at the following to get the number of repetitions automatically!\n"
            m1 <- paste(m1, paste(where, collapse="\n"), collapse="")
            m1 <- paste(m1, "\n The existing files in this folder:\n")
            m1 <- paste(m1, paste(list.files(dir), collapse="\n"))
            stop(m1)
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
    return(result)
}
