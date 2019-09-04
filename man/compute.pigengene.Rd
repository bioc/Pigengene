\name{compute.pigengene}
\alias{compute.pigengene}

\title{
Computes the eigengenes}

\description{
This function takes as input the expression data and module assignments, and 
computes an eigengene for each module using PCA.}

\usage{
compute.pigengene(Data, Labels, modules, saveFile = "pigengene.RData", 
    selectedModules = "All", amplification = 5, doPlot = TRUE,
    verbose = 0, dOrderByW = TRUE, naTolerance=0.05)
}


\arguments{
  \item{Data}{
  A matrix or data frame containing the training expression data, with genes 
  corresponding to columns and rows corresponding to samples. Rows and columns 
  must be named.}

  \item{Labels}{
  A (preferably named) vector containing the Labels (condition types) for the training Data.
  Names must agree with rows of \code{Data}.}

\item{modules}{A numeric vector, named by \code{genes}, 
    that reports the module (clustering) assignments.}

  \item{saveFile}{
   The file to save the results. \code{NULL} will disable saving, and
   thus requires \code{doPlot} to be \code{FALSE}.}

  \item{selectedModules}{
  A numeric vector determining which modules to use, or set to 
  "All" (default) to include every module.}

  \item{amplification}{
  An integer that controls the number of repeats for each condition.
  The number of all samples roughly will be multiplied by this factor 
  after oversampling. See \code{\link{balance}}.}

  \item{doPlot}{
  Boolean determining whether heatmaps of expression of eigengenes 
  should be ploted and saved.}

  \item{verbose}{
  The integer level of verbosity. 0 means silent and higher values    
  produce more details of computation.}

 \item{dOrderByW}{
  If \code{TRUE}, the genes will be ordered in the csv file based on their
  absolute weight in the corresponding module.}

  \item{naTolerance}{
   Upper threshold on the fraction of entries per gene that
   can be missing. Genes with a larger fraction of missing
   entries are ignored. For genes with smaller fraction of NA
   entries, the missing values are imputed from their average
   expression in the other samples.  
   See \code{\link{check.pigengene.input}}.}

}

\details{
Rows of \code{Data} are oversampled using \code{\link{balance}} so that 
each condition has roughly the same number of samples. 
\code{\link[WGCNA]{moduleEigengenes}} computes an eigengene for each module
using PCA.
}

\value{
An object of \code{\link{pigengene-class}}.}

\references{
  Large-scale gene network analysis reveals the significance of extracellular 
  matrix pathway and homeobox genes in acute myeloid leukemia,
  Foroushani A, Agrahari R, Docking R, Karsan A, and Zare H.
  In preparation.}

\author{Habil Zare and Amir Foroushani}

\seealso{
\code{\link{Pigengene-package}},
\code{\link{one.step.pigengene}}, \code{\link{wgcna.one.step}},
\code{\link{make.decision.tree}}, \code{\link[WGCNA]{moduleEigengenes}}}

\examples{
## Data:
data(aml)
data(mds)
data(eigengenes33)
d1 <- rbind(aml,mds)
Labels <- c(rep("AML",nrow(aml)),rep("MDS",nrow(mds)))
names(Labels) <- rownames(d1)
modules33 <- eigengenes33$modules[colnames(d1)]
## Computing:
pigengene <- compute.pigengene(Data=d1, Labels=Labels, modules=modules33,
   saveFile="pigengene.RData", doPlot=TRUE, verbose=3)
class(pigengene)
plot(pigengene, fontsize=12)
}

\keyword{ cluster }
