\name{Pigengene-package}
\alias{Pigengene-package}
\docType{package}
\title{Infers robust biological signatures from gene expression data}

\description{
Pigengene identifies gene modules (clusters), computes an eigengene for each 
module, and uses these biological signatures as features for classification.
The resulting biological signatures are very robust with respect to the 
profiling platform. For instance, if Pigenegene computes a biological signature
using a microarray dataset, it can infer the same signature in an RNA Seq
dataset such that it is directly comparable across the two datasets.}

\details{
\tabular{ll}{
Package: \tab Pigengene\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2016-04-25\cr
License: \tab GPL (>= 2)\cr
}

The main function is \code{\link{one.step.pigengene}} which requires a gene 
expression profile and the corresponding conditions (types).
Individual functions are provided to facilitate running the pipeline in a 
customized way. Also, the inferred biological signatures (computed eigengenes)
are useful for other supervised or unsupervised analyses.

In most functions of this package, eigenegenes are computed or used as robust
biological signatures. Briefly, each eigengene is a weighted average of the
expression of all genes in a module (cluster), where the weights are adjusted 
in a way that the explained variance is maximized. 
}

\author{
Amir Foroushani, Habil Zare, and Rupesh Agrahari

Maintainer: Habil Zare <zare@txstate.edu>
}

\references{
  Large-scale gene network analysis reveals the significance of extracellular 
  matrix pathway and homeobox genes in acute myeloid leukemia,
  Foroushani A, Agrahari R, Docking R, Karsan A, and Zare H.
  In preparation.}

\keyword{ package }
\keyword{ datasets }
\keyword{ documentation }
\keyword{ models }

\seealso{
\code{\link[Pigengene]{Pigengene-package}},
\code{\link{one.step.pigengene}}, \code{\link{compute.pigengene}}, 
\code{\link[WGCNA]{WGCNA-package}}}


\examples{
data(aml)
data(mds)
d1 <- rbind(aml,mds)
Labels <- c(rep("AML",nrow(aml)),rep("MDS",nrow(mds)))
names(Labels) <- rownames(d1)
p1 <- one.step.pigengene(Data=d1,saveDir='pigengene', bnNum=10, verbose=1,
   seed=1, Labels=Labels, toCompact=FALSE, doHeat=FALSE)
plot(p1$c5treeRes$c5Trees[["34"]])
## See pigengene for results.
}
