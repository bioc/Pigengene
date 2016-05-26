\name{get.fitted.leaf}
\alias{get.fitted.leaf}

\title{Returs the leaf for each sample}

\description{
Taking as input a tree and data, this function determines the leaf
each sample will fall in.}

\usage{
get.fitted.leaf(c5Tree, inpDTemp, epsi = 10^(-7))
}

\arguments{
  \item{c5Tree}{
    A decision tree of class \code{C50} that uses module eigengenes, 
    or \code{NULL}.  If \code{NULL}, expression plots for all modules are
    created.}

  \item{inpDTemp}{
    The possibly new data matrix with samples on rows}

  \item{epsi}{
    A small perturbation to resolve the boundary issue}
}

\value{
A numeric vector of node indices named by samples (rows of \code{inpDTemp})}

\author{Amir Foroushani}

\note{
This function is tricky because C50 uses a global variable.}

\seealso{
\code{\link{Pigengene-package}},
\code{\link{make.decision.tree}}, \code{\link{compact.tree}}, 
\code{\link{compute.pigengene}}, \code{\link{module.heatmap}},
\code{\link{get.used.features}}, \code{\link{preds.at}}
}

\examples{
     ## Data:
     data(aml)
     data(mds)
     data(pigengene)
     d1 <- rbind(aml,mds)

     ## Fiting the trees:
     trees <- make.decision.tree(pigengene=pigengene, Data=d1,
     saveDir="trees", minPerLeaf=15, doHeat=FALSE,verbose=3,
       toCompact=FALSE)
     f1 <- get.fitted.leaf(c5Tree=trees$c5Trees[["15"]], 
       inpDTemp=pigengene$eigengenes)
}

\keyword{ tree }
\keyword{ utilities }
