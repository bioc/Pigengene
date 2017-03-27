\name{module.heatmap}
\alias{module.heatmap}

\title{Plots heatmaps for modules}

\description{
This function takes as input a tree and an object from 
\code{\link{pigengene-class}} and per any module used in the tree,
it plots one gene expression heatmap.}

\usage{
module.heatmap(c5Tree, pigengene, saveDir, testD = NULL,
  testL = NULL, pos = 0, verbose=0, doAddEigengene=TRUE, ...)
}

\arguments{
  \item{c5Tree}{
    A decision tree of class \code{C50} that uses module eigengenes, 
    or \code{NULL}. If \code{NULL}, expression plots for
    all modules are created.}

  \item{pigengene}{
    A object of \code{\link{pigengene-class}}, output of 
    \code{\link{compute.pigengene}}}

  \item{saveDir}{
    Directory to save the plots}

  \item{testD, testL}{
    Optional. The matrix of (independent) test
    expression data and the corresponding vector of labels}

  \item{pos}{
    Number of genes to discard. Interpreted the same way as in
    \code{\link{compact.tree}} ad \code{\link{preds.at}}}

  \item{verbose}{
    The integer level of verbosity. 0 means silent and higher values    
    produce more details of computation.}

  \item{doAddEigengene}{
    If \code{TRUE}, the eigengene of each module will be added to the
    corresponding heatmap.}
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
    
     ## Fiting the trees:
     trees <- make.decision.tree(pigengene=pigengene, Data=d1,
       saveDir="trees", minPerLeaf=14:15, doHeat=FALSE,verbose=3,
       toCompact=15)

     ## Plotting:
     module.heatmap(c5Tree=trees$c5Trees[["15"]], pigengene=pigengene, 
         saveDir="heatmaps", pos=0, verbose=1)

}

\keyword{ cluster }
\keyword{ tree }
\keyword{ utilities }
