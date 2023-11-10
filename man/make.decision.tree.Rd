\name{make.decision.tree}
\alias{make.decision.tree}

\title{Creates a decision tree to classify samples using the eigengenes values}

\description{
A decision tree in \code{\link{Pigengene-package}} uses module eigengenes to
build a classifier that distinguishes the different classes.
Briefly, each eigengene is a weighted average of the expression of
all genes in the module, where the weight of each gene corresponds
to its membership in the module.}

\usage{
make.decision.tree(pigengene, Data, 
  Labels = structure(pigengene$annotation[rownames(pigengene$eigengenes),
          1], names = rownames(pigengene$eigengenes)),
  testD = NULL, testL = NULL, selectedFeatures = NULL,
  saveDir = "C5Trees", minPerLeaf = NULL, useMod0 = FALSE, 
  costRatio = 1, toCompact = NULL, noise = 0, noiseRepNum = 10, doHeat=TRUE,
  verbose = 0, naTolerance=0.05)
}

\arguments{
  \item{pigengene}{
    The pigengene object that is used to build the decision tree. 
    See \code{\link{pigengene-class}}.}

  \item{Data}{
    The training expression data}

  \item{Labels}{
    Labels (condition types) for the (training) expression data. It is a
    named vector of characters. \code{Data} and \code{pigengene} will be
    subset according to these names.}

  \item{testD}{
    The test expression data, for example, from an independent dataset. 
    Optional.}

  \item{testL}{
     Labels (condition types) for the (test) expression data. Optional.}

  \item{selectedFeatures}{
    A numeric vector determining the subset of eigengenes
    that should be used as potential predictors. By default
    ("All"), eigengenes for all modules are considered.  
    See also \code{useMod0}.}

  \item{saveDir}{
    Where to save the plots of the tree(s).}

  \item{minPerLeaf}{
    Vector of integers. For each value, a tree will be
    built requiring at least that many nodes on each leaf. By
    default (\code{NULL}), several trees are built, one for each
    possible value between 2 and 10 percent of the number of samples.}

  \item{useMod0}{
    Boolean. Wether to allow the tree(s) to use the eigengene of
    module 0, which corresponds to the set of outlier, as a proper predictor.}

  \item{costRatio}{
    A numeric value effective only for 2 groups classification.
    The default value (1) considers the misclassification of both
    conditions as equally disadvatageous. Change this value to a
    larger or smaller value if you are more interested in the
    specificity of predictions for condition 1 or condition 2, respectively.}

  \item{toCompact}{
    An integer. The tree with this \code{minPerLeaf} value will
    be compacted (shrunk). Compacting in this context means
    reducing the number of required genes for the calculation of
    the relevant eigengenes and making the predictions using the
    tree. If \code{NULL} (default), the (persumably) most general
    proper tree (corresponding to the largest value in the
    \code{minPerLeaf} vector for which a tree could be
    constructed) is compacted. Set to \code{FALSE} to turn off compacting.}

  \item{noise, noiseRepNum}{
    For development purposes only. These parameters
    allow investigating the effect of gaussian noise in the
    expression data on the accurracy of the tree for test data.}

  \item{doHeat}{
    Boolean. Set to \code{FALSE} not to plot the heatmaps 
    for faster comoutation.}
  
  \item{verbose}{
    The integer level of verbosity. 0 means silent and higher values    
    produce more details of computation.}

  \item{naTolerance}{
    Upper threshold on the fraction of entries per gene that
    can be missing. Genes with a larger fraction of missing
    entries are ignored. For genes with smaller fraction of NA
    entries, the missing values are imputed from their average
    expression in the other samples.  
    See \code{\link{check.pigengene.input}}.}
}

\details{
This function passes the inut eigengenes and appropriate arguments
to \code{\link[C50]{C5.0}} function from \code{C50} package.

The effect of test data:
Only when both \code{testD} and \code{testL} are provided,
the test data will be used for a) compacting the trees, 
b) plotting heatmaps of expression of genes in the
compacted and full trees, and
c) the noise analysis. 
If either of \code{testD} or \code{testL} is \code{NULL}, then
\code{Data} and \code{Labels} are instead used for these purposes.
}

\value{
A list with following elements:

\item{call}{
  The call that created the results}
\item{c5Trees}{
  A list, with one element of class \code{C5.0} for each attempted
  \code{minNodesperleaf} value.  The list is named with the
  corresponding values as characters. An extra \code{info} element is
  added that includes information on the performance of the tree.}
\item{minPerLeaf}{
  A numeric vector enumerating all of the attempted
  minPerLeaf values.}
\item{compacted}{
  The full output of \code{\link{compact.tree}} function if 
  \code{toCompact} is not \code{FALSE}}
\item{heat}{
  The output of  \code{\link{module.heatmap}} function for the full tree if
  \code{doHeat} is not \code{FALSE}}
\item{heatCompact}{
  The output of  \code{\link{module.heatmap}} function for the compacted
  tree if \code{toCompact} is not \code{FALSE}}
\item{noisy}{
  The full output of  \code{noise.analysiy} function if \code{noise} is
  not 0. For development and evaluation purposes only.}
\item{leafLocs}{
  A matrix reporting the leaf for each sample on 1 row.  The columns are named
  according to the correspoding \code{minNodesperleaf} value.}
\item{toCompact}{
  Echos the \code{toCompact} input argument}
\item{costs}{
  The cost matrix}
\item{saveDir}{
  The directory where plots are saved in}
}

\note{
For faster computation in an initial, explanatory run, turn off
compacting, which can take a few minutes, with \code{toCompact=FALSE}.}

\seealso{
\code{\link{Pigengene-package}}, \code{\link{compute.pigengene}},
\code{\link{compact.tree}}, \code{\link[C50]{C5.0}},
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
       toCompact=15)
}

\keyword{ tree }
\keyword{ classif }
