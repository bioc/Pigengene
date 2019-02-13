get.sorted.members <- function(modID, pigengene, modules ){
    numModID <- as.numeric(gsub("ME", "", modID))
    ## If there was no "ME" in modID, the old version would look up the wrong column 
    charModID <- paste("ME", numModID, sep='')
    Q1 <- pigengene$membership[match(names(modules)[which(modules==numModID)], rownames(pigengene$membership)), charModID]
    return(Q1[order(rank(abs(Q1)))])
}





