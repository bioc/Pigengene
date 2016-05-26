\name{get.genes}
\alias{get.genes}
\title{List the (most relevant) genes for a decision tree.}

\description{
This function returns all genes that are left after shrinking (compacting ) 
a given tree. If \code{enhance} is set to \code{TRUE}, it makes sure that 
the output contains at least two genes from each used module.}

\usage{
get.genes(c5Tree = NULL, pigengene = NULL, queue = NULL, modules = NULL, pos=0,
  enhance = TRUE)
}

\arguments{
  \item{queue}{
  A character vector. The membership queue for a decsision tree.}

  \item{pos}{
  Number of genes that are considered from removal. Same interpretation as in
  \code{\link{preds.at}}}

  \item{enhance}{
  If \code{enhance} is set to \code{TRUE}, the function makes sure that 
  the output contains at least two genes from each used module. 
  Otherwise, exactly the \code{pos} first elements of the queue are removed
  from consideration.}

  \item{modules}{
  Named character vector listing the module assignments.}

  \item{c5Tree}{
  A decision tree of class \code{C50}.}

  \item{pigengene}{
  An object in \code{\link{pigengene-class}}, usually created by 
  \code{\link{compute.pigengene}}.}
}

\details{
This function needs \code{modules} and \code{queue}, or alternatively,  
\code{c5Tree} and \code{pigengene}.}

\value{
  A character vector containing the names of the genes involved in the
  modules whose eigengenes are used in the tree. If \code{pos > 0}, 
  the first \code{pos} such genes with lowest absolute membership in their 
  respective modules are filtered.
}

\seealso{
  \code{\link{Pigengene-package}},
  \code{\link{compact.tree}},\code{\link{preds.at}},
  \code{\link{get.used.features}}, \code{\link{make.decision.tree}}
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
     g1 <- get.genes(c5Tree=trees$c5Trees[["15"]],pigengene=pigengene)
}

\keyword{ tree }
