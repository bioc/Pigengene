\name{compact.tree}
\alias{compact.tree}

\title{Reduces the number of genes in a decision tree}

\description{
In a greedy way, this function removes the genes with smaller
weight one-by-one, while assessing the accuracy of the predictions
of the resulting trees.}

\usage{
compact.tree(c5Tree, pigengene, Data=pigengene$Data, Labels=pigengene$Labels,
  testD=NULL, testL=NULL, saveDir=".", verbose=0)
}

\arguments{
  \item{c5Tree}{
    A decision tree of class \code{C50} that uses module eigengenes, 
    or \code{NULL}. If \code{NULL}, If \code{NULL}, expression plots for
    all modules are created.}

  \item{pigengene}{
    A object of \code{\link{pigengene-class}}, output of 
    \code{\link{compute.pigengene}}}

  \item{Data}{
    A matrix or data frame containing the expression data, with
    genes corresponding to columns and rows corresponding to
    samples. Rows and columns must be named.}

  \item{Labels}{
    Labels (condition types) for the (training) expression data.  It is a
    named vector of characters. \code{Data}  will be subset according to
    these names.}

  \item{testD}{
    The test expression data, for example, from an independent 
    dataset. Optional.}

  \item{testL}{
    Labels (condition types) for the (test) expression data. Optional.}

  \item{saveDir}{
    Where to save the plots of the tree(s)}

  \item{verbose}{
  Integer level of verbosity. 0 means silent and higher values    
  produce more details of computation.} 

}

\value{
A list with following elements is invisibly returned:
\item{call}{
  The call that created the results}
\item{predTrain}{
  Prediction using projected data without compacting}
\item{predTrainCompact}{
  Prediction after compacting}
\item{genes}{
  A character vector of all genes in the full tree before compacting}
\item{genesCompacted}{
  A character vector of all genes in the compacted tree}
\item{trainErrors}{
  A matrix reporting errors on the train data.  The rows are
  named according to the number of removed genes. Each column
  reports the number of misclassified samples in one condition
  (type) except the last column that reports the total.}
\item{testErrors}{
  A matrix reporting errors on the test data similar to \code{trainErrors}}
\item{queue}{
  A numeric vector named by all genes contributing to the full
  tree before compacting. The numeric values are weights
  increasingly ordered by absolute value.}
\item{pos}{
  The number of removed genes}
\item{txtFile}{
  Confusion matrices and other details on compacting are reported 
  in this text file}
}

\references{
     Large-scale gene network analysis reveals the significance of
     extracellular matrix pathway and homeobox genes in acute myeloid
     leukemia, Foroushani A, Agrahari R, Docking R, Karsan A, and Zare
     H.  In preparation.

     Gene shaving as a method for identifying distinct sets of genes
     with similar expression patterns, Hastie, Trevor, et al. Genome
     Biol 1.2 (2000): 1-0003.
}

\seealso{
\code{\link{Pigengene-package}}, \code{\link{compute.pigengene}},
\code{\link{make.decision.tree}}, \code{\link[C50]{C5.0}},
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
     saveDir="trees", minPerLeaf=14:15, doHeat=FALSE,verbose=3,
     toCompact=FALSE)
     c1 <- compact.tree(c5Tree=trees$c5Trees[["15"]], pigengene=pigengene,
         saveDir="compacted", verbose=1)

}

\keyword{ tree }
\keyword{ classif }
