draw.module.sizes <- function(
    net=NULL, plotFile=NULL, verbose=1) 
{
    outNum <- length(net$colors[net$colors==0])
    mainText <- paste("Size of ", length(unique(net$colors)) - 1, 
                      " modules with", outNum, " outliers", sep="")
    if(!is.null(plotFile)){
        dir.create(path=dirname(plotFile), recursive=TRUE, showWarnings=FALSE)
        png(filename=plotFile)
    }
    barplot(table(net$colors[net$colors > 0]), col="blue", 
            main=mainText, ylab="Number of genes", xlab="Module index", 
            cex.lab=1.3)
    if(!is.null(plotFile)){
        dev.off()
    }
}
