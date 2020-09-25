\name{module.heatmap}
\alias{module.heatmap}

\title{Plots heatmaps for modules}

\description{
This function takes as input a tree and an object from 
\code{\link{pigengene-class}} and per any module used in the tree,
it plots one gene expression heatmap. Alternatively, it can plot
a heatmap for every module in the given \code{pigengene} object.}

\usage{
module.heatmap(c5Tree=NULL, pigengene, Data=NULL, mes=NULL, saveDir, testD = NULL,
  testL = NULL, pos = 0, verbose=0, doAddEigengene=TRUE, scalePngs=1, ...)
}

\arguments{
  \item{c5Tree}{
    A decision tree of class \code{C50} that uses module eigengenes, 
    or \code{NULL}. If \code{NULL}, expression plots for all modules are
    created.}

  \item{pigengene}{
    A object of \code{\link{pigengene-class}}, output of 
    \code{\link{compute.pigengene}}}

  \item{Data}{
    The gene expression matrix with genes on columns and samples on
    rows. Set it to NULL to use the data in the \code{pigengene}
    input argument.}

  \item{mes}{
    A character vector that determines which modules to plot, e.g.,
    c("ME3","ME5"). Set it to \code{NULL} to plot a heatmap for every
    module. This argument will be ignored if \code{c5Tree} is not \code{NULL}.}

  \item{saveDir}{
    Directory to save the plots}

  \item{testD, testL}{
    Optional. The matrix of (independent) test
    expression data and the corresponding vector of labels}

  \item{pos}{
    Number of genes to discard. Interpreted the same way as in
    \code{\link{compact.tree}} and \code{\link{preds.at}}}

  \item{verbose}{
    The integer level of verbosity. 0 means silent and higher values    
    produce more details of computation.}

  \item{doAddEigengene}{
    If \code{TRUE}, the eigengene of each module will be added to the
    corresponding heatmap.}

  \item{scalePngs}{
    If not 1, the size of pngs will be adjusted using this
    parameter. A typical value would be 7.}
  
  \item{\dots}{
    Additional arguments. Passed to \code{\link{pheatmap.type}}}

}

\value{
A list of:
  \item{call}{
    The call that created the results}
  \item{saveDir}{
    An echo of the input argument determining 
    where the plots are saved}
}

\seealso{
\code{\link{Pigengene-package}},
\code{\link{make.decision.tree}}, \code{\link{compact.tree}}, 
\code{\link{compute.pigengene}}
}

\examples{
     ## Data:
     data(aml)
     data(mds)
     data(pigengene)
     d1 <- rbind(aml,mds)

     ## Plotting the heatmaps of all modules:
     module.heatmap(pigengene=pigengene, saveDir="heatmaps", pos=0, verbose=1)
    
     ## Fiting the trees:
     trees <- make.decision.tree(pigengene=pigengene, Data=d1,
       saveDir="trees", minPerLeaf=14:15, doHeat=FALSE,verbose=3,
       toCompact=15)

     ## Plotting the heatmaps of only the modules in the tree:
     module.heatmap(c5Tree=trees$c5Trees[["15"]], pigengene=pigengene, 
         saveDir="treeHeatmaps", pos=0, verbose=1)

}

\keyword{ cluster }
\keyword{ tree }
\keyword{ utilities }
