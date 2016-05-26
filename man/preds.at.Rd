\name{preds.at}
\alias{preds.at}

\title{Prediction using a possibly compacted tree}

\description{
A decision tree in Pigengene uses module eigengenes to build a classifier that
distincuishes two or more classes. Each eigengene is a weighted average of the
expression of all genes in the module, where the weight of each gene corresponds
to its membership in the module. Each modules might contain dozens to hundreds
of genes, and hence the final classifier might depend on the expression of 
a large number of genes. In practice, it can be desireable to reduce the number
of necessary genes used by a decision tree. This function is helpful in observing
changes to the classification output after removing genes with lower weights
membership. It determines how a given decision tree would classify the expression
data after removing a certain number of genes from consideration.}

\usage{
preds.at(c5Tree, pigengene, pos=0, Data)
}

\arguments{
  \item{c5Tree}{
    A decision tree that uses eigengenes from the pigengene
    object to classify the samples from the expression data.}

  \item{pigengene}{
    A object of \code{\link{pigengene-class}}, output of 
    \code{\link{compute.pigengene}}}

  \item{pos}{
    Number of genes to be removed from the consideration. Genes
    are removed in ascending order of their absolute weight in
    the relevant modules.  If 0 (default), the prediction will be
    done without compacting.}

  \item{Data}{
    The expression possibly new data used for classification}
}

\value{
A list with following components:

\item{predictions}{
  The vector of predictions after neglecting \code{pos} number of genes}
\item{eigengenes}{
  The values for the eigenges after neglecting \code{pos} number of genes}
}

\seealso{
\code{\link{Pigengene-package}}, \code{\link{pigengene-class}},
\code{\link{make.decision.tree}}, \code{\link{compact.tree}}, 
\code{\link{compute.pigengene}}, \code{\link{module.heatmap}},
\code{\link{get.used.features}}, \code{\link{get.fitted.leaf}},
\code{\link{Pigengene-package}}
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
     preds1 <- preds.at(c5Tree=trees$c5Trees[["15"]], pigengene=pigengene, 
       pos=0, Data=d1)
}

\keyword{ tree }
