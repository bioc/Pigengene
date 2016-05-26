unique.gene.names <- function (ids, inputType="ENTREZIDat"){
    result <- c()
    if (length(ids)==0) 
        return(ids)
    geneRes <- gene.mapping(ids, inputType=inputType)
    geneVector <- geneRes[, "output2"]
    geneIds <- geneRes[, "input"]
    resGeneVector <- c(rep(x=NA, length(geneIds)))
    for (geneI in 1:length(geneVector)) {
        indSameGene <- which(geneVector==geneVector[geneI])
        if (length(indSameGene) > 1) {
            for (indices in indSameGene) {
                resGeneVector[indices] <- paste(geneVector[indices], 
                                                geneIds[indices], sep="_")
            }
        }
    }
    naInd <- which(is.na(resGeneVector))
    for (kIndex in naInd) {
        resGeneVector[kIndex] <- geneVector[kIndex]
    }
    naInd <- which(is.na(resGeneVector))
    for (kIndex in naInd) {
        resGeneVector[kIndex] <- geneIds[kIndex]
    }
    names(resGeneVector) <- ids
    result[["geneVector"]] <- resGeneVector
    result[["ids"]] <- ids
    result[["geneIds"]] <- geneIds
    return(result)
}
