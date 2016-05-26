get.fitted.leaf <- function(c5Tree, inpDTemp, epsi=10^(-7)){
    ## temporarily overwrite the global inpD, get the leafs on inpDTemp
    ##^as input, restore original inpD
    ## epsi: epsilon, set to 1E-7 to handle cases on the boundary.
    if(class(c5Tree) !="C5.0")
        stop("The class of c5Tree argument must be 'C5.0' !")
    inpIGlobalName <- as.character(c5Tree$call)[3]
    inpDbk <- get(inpIGlobalName) ## inpD
    if(class(inpDbk) != "data.frame")
        stop(paste("Could not get the global variable ", inpIGlobalName,"!"))
    if(! "Labels" %in% colnames(inpDTemp))
        inpDTemp <- cbind(inpDTemp,Labels="Labels")
    if(any(! colnames(inpDbk) %in% colnames(inpDTemp))){
        stop(paste("Colnames of inpDTemp do not match witht the",
                   inpIGlobalName,"!"))
    }
    if(epsi>0){
        notypeCols <- grep("ME", colnames(inpDTemp), value=TRUE)
        inpDTemp[, notypeCols] <- inpDTemp[, notypeCols]+epsi
    }
    assign('inpD', inpDTemp , envir=parent.env(parent.frame()), inherits=TRUE)
    qq <- eval(parse(text=(as.character(partykit::as.party(c5Tree))[3])))[[1]]
    assign('inpD', inpDbk , envir=parent.env(parent.frame()), inherits=TRUE)
    names(qq) <- rownames(inpDTemp)
    return(qq)
}
