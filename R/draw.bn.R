draw.bn <- function(
    BN, plotFile=NULL, inputType="ENTREZIDat", edgeColor="blue", 
    DiseaseCol="darkgreen", DiseaseFill="red", 
    DiseaseChildFill="pink", ##"olivedrab1"
    nodeCol="darkgreen", nodeFill="yellow", moduleNamesFile=NULL, 
    mainText=NULL, nodeFontSize=14 * 1.1, verbose=0)
{
    ## BN: an object of class BN.
    ## plotFile: set to NULL not to save the plot.
    ## nodeFontSize: the default is 14 in Rgraphviz package.
    message.if(me="Drawing BN ...", verbose=verbose)
    result <- list()
    result[["call"]] <- match.call()
    a0 <- Rgraphviz::getDefaultAttrs()
    result[["BN"]] <- BN
    if (!is.null(moduleNamesFile)) {
        BN <- rename.nodes(BN=BN, inputType=inputType, moduleNamesFile=moduleNamesFile)$bn
        result[["renamedBN"]] <- BN
    }
    gr <- graphviz.plot(BN)
    nodeNames <- graph::nodes(gr)
    edgesList <- graph::edges(gr)
    DiseaseChildren <- edgesList[["Disease"]] ## A character vector.
    nNodes <- length(nodeNames)
    generalNodesAttrList <- list(shape="ellipse", fixedsize=TRUE)
    generalEdgesAttrList <- list(color=edgeColor)
    nAttrs <- list()
    eAttrs <- list()
    attrs <- list(node=generalNodesAttrList, edge=generalEdgesAttrList)
    nAttrs$color <- rep(nodeCol, nNodes)
    nAttrs$fillcolor <- rep(nodeFill, nNodes)
    nAttrs$fillcolor[which(nodeNames %in% DiseaseChildren)] <- DiseaseChildFill
    nAttrs$fontsize <- rep(nodeFontSize, nNodes)
    nAttrs$width <- rep(as.numeric(a0$node$width) * (nodeFontSize/14)^9, 
                        nNodes)
    nAttrs$height <- rep(as.numeric(a0$node$height) * (nodeFontSize/14)^8, 
                         nNodes)
    nAttrs <- lapply(nAttrs, function(x) {
        names(x) <- nodeNames
        x
    })
    nAttrs$fillcolor["Disease"] <- DiseaseFill
    nAttrs$color["Disease"] <- DiseaseCol
    nAttrs$fontsize["Disease"] <- nodeFontSize * 1.5
    if (!is.null(plotFile)) 
        png(plotFile, width=480 * 3, height=480 * 3)
    plotted <- graph::plot(gr, attrs=attrs, nodeAttrs=nAttrs, edgeAttrs=eAttrs, 
                           main=mainText)
    if (!is.null(plotFile)) {
        dev.off()
        message.if(paste("The BN graph was plotted in:", plotFile), verbose=verbose-1)
    }
    result[["gr"]] <- gr
    result[["plotted"]] <- plotted
    result[["plotFile"]] <- plotFile
    return(result)
}
