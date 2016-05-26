rename.nodes <- function(
    BN=NULL, bnetIndv=NULL, saveFile=NULL, 
    inputType="ENTREZIDat", 
    moduleNamesFile=NULL, colName="MName")
{
    ## colName: The name of the column in moduleNamesFile that has the new module names.
    result <- list()
    origNodeList <- c()
    mappingTable <- NULL
    if (!is.null(moduleNamesFile)) {
        mappingTable <- as.matrix(read.csv(moduleNamesFile))
        mapping <- mappingTable[, colName]
        names(mapping) <- mappingTable[, "ModuleID"]
    }
    get.new.nodeList <- function(origNodeList) {
        if (!is.null(mappingTable)) {
            nodeList <- mapping[origNodeList]
        }
        else {
            nodeList <- unique.gene.names(ids=origNodeList, 
                                          inputType=inputType)$geneVector
        }
        return(nodeList)
    }
    if (is.null(BN)) {
        if (is.null(bnetIndv)) {
            stop("Both BN and bnet are null!")
        }
        else {
            arcsBN <- bnetIndv
            origNodeList <- unique(as.character(arcsBN[, c("from", 
                                                           "to")]))
            nodeList <- get.new.nodeList(origNodeList)
            bnObj <- empty.graph(nodes=nodeList)
        }
    }
    else {
        bnObj <- BN
        arcsBN <- bnlearn::arcs(BN)
        origNodeList <- bnlearn::nodes(BN)
        nodeList <- get.new.nodeList(origNodeList)
    }
    bnlearn::nodes(bnObj) <- nodeList
    bnlearn::arcs(bnObj) <- cbind(nodeList[arcsBN[, "from"]], nodeList[arcsBN[, "to"]])
    print("BN object created with new node names.")
    if (length(origNodeList) !=length(nodes(bnObj))) 
        stop("Number of nodes should not change by renaming!")
    result[["inputBN"]] <- BN
    result[["bnetIndv"]] <- bnetIndv
    result[["bn"]] <- bnObj
    result[["mappingTable"]] <- mappingTable
    renamed <- result
    if (!is.null(saveFile)) {
        save.if(renamed, file=saveFile)
        print(paste("renamed was saved in:", saveFile))
    }
    return(renamed)
}
