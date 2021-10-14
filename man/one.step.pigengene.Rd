\name{one.step.pigengene}
\alias{one.step.pigengene}

\title{
Runs the entire Pigengene pipeline
}

\description{
Runs the entire Pigengene pipeline, from gene expression to compact decision 
trees in a single function. It identifies the gene modules using coexpression
network analysis, computes eigengenes, learns a Bayesian network, fits 
decision trees, and compact them.}

\usage{
one.step.pigengene(Data, saveDir="Pigengene", Labels, testD=NULL, 
  testLabels=NULL, doBalance=TRUE, RsquaredCut=0.8, costRatio=1,
  toCompact=FALSE, bnNum=0, bnArgs=NULL, useMod0=FALSE, mit="All",
  verbose=0, doHeat=TRUE, seed=NULL, dOrderByW=TRUE, naTolerance=0.05, 
  idType="ENTREZID", pathwayDb=NULL, OrgDb=org.Hs.eg.db)
}

\arguments{
  \item{Data}{
  A matrix or data frame (or list of matrices or data frames) containing
  the training expression data, with genes  corresponding to columns and
  rows corresponding to samples. Rows and columns 
  must be named. For example, from RNA-Seq data, log(RPKM+1) can be used.}
  
  \item{Labels}{
  A (preferably named) vector containing the Labels (condition types)
  for the training Data. Or, if Data is a list, a list of label vectors
  corresponding to the data sets in Data.
  Names must agree with rows of \code{Data}.}

  \item{saveDir}{
  Directory to save the results.}

  \item{testD}{
  Test expression data with syntax similar to \code{Data},
  possibly with different rows and columns. This argument is optional
  and can be set to \code{NULL} if test data are not available.}

  \item{testLabels}{
  A (preferably named) vector containing the Labels (condition types) for
  the test Data. This argument is optional and can be set to \code{NULL}
  if test data are not available.}

  \item{doBalance}{
  Boolean. Whether the data should be oversampled before identifying the modules
  so that each condition contribute roughly the same number of samples
  to clustering.}

  \item{RsquaredCut}{
  A threshold in the range [0,1] used to estimate the power. A higher value
  can increase power. For technical use only. See \code{\link[WGCNA]{pickSoftThreshold}}
  for more details. A larger value generally leads to more modules.}

  \item{costRatio}{
  A numeric value, the relative cost of misclassifying a sample from the first
  condition vs. misclassifying a sample from the second condition.}
  
  \item{toCompact}{
  An integer value determining which decision tree to shrink. 
  It is the minimum number of genes per leaf imposed when fitting the tree.
  Set to \code{FALSE} to skip compacting, to \code{NULL} to automatically 
  select the maximum value.}

  \item{bnNum}{
  Desired number of bootstraped Baysian networks. 
  Set to \code{0} to skip BN learning.}

  \item{bnArgs}{
  A list of arguments passed to \code{\link{learn.bn}} function.}

  \item{useMod0}{
  Boolean, whether to allow module zero (the set of outliers) to be used
  as a predictor in the decision tree(s).}

  \item{mit}{
  The "module identification type", a character vector determining the reference
  conditions for clustering. If 'All' (default), clustering is performed using 
  the entire data regardless of condition.}

  \item{verbose}{
  The integer level of verbosity. 0 means silent and higher values    
  produce more details of computation.} 

  \item{doHeat}{
  If \code{TRUE} the heatmap of expression of genes in the modules that 
  contribute to the the tree will be plotted.}

  \item{seed}{
  Random seed to ensure reproducibility.}

  \item{dOrderByW}{
  If \code{TRUE}, the genes will be ordered in the csv file based on their
  absolute weight in the corresponding module.}
  
  \item{naTolerance}{
  Upper threshold on the fraction of entries per gene that
  can be missing. Genes with a larger fraction of missing
  entries are ignored. For genes with smaller fraction of NA
  entries, the missing values are imputed from their average
  expression in the other samples. See \code{\link{check.pigengene.input}}.}

  \item{idType}{
  A string describing the type of input gene ID e.g., "ENTREZID", "REFSEQ", 
  "SYMBOL".}

  \item{pathwayDb}{
  A character vector determining which enrichment database to be used by the
  \code{\link{get.enriched.pw}} function e.g., "GO", "KEGG", "REACTOME", or
  "NCG". Set to \code{NULL} to skip the pathway enrichment analysis.}

  \item{OrgDb}{
  The reference data base to be used. Use e.g. \code{org.Ce.eg.db} for 'Celegans'
  when analysing Celegans data. If \code{OrgDb} is not \code{NULL}, \code{Org} 
  must be \code{NULL}.}
}

\details{
This is the main function of the package Pigengene and performs several
steps: First, modules are identified in the training expression data,
according to mit argument i.e.  based on coexpression behaviour
in the corresponding conditions. Set it to "All" to use all training data for 
this step regardless of the condition. Then, if a list of data frames is
provided in Data, similarity networks on the data sets are computed and
combined into one similarity network for the union of nodes across data
sets. Then, the eigengenes for each module 
and each sample are calculated, where the expression of an eigengene of a
module in a sample is the weighted average of the expression of the genes in 
that module in the sample. Technically, an eigengene is the first principal 
component of the gene expression in a module. PCA ensures that the maximum 
variance accross all the training samples is explained by the eigengene. 
Next, (optionally --if bnNum is set to a value greater than 0),
several bootstrapped Bayesian networks are learned and combined into a
consensus network, in order to detect and illustrate the 
probabilistic dependencies between the eigengenes and the disease subtype. 
Next, decisision tree(s) are built that use the module eigengenes, or
a subset of them, to  distinguish the classes (\code{Labels}). 
The accurracy of trees is assessed on the train and (if provided) test data. 
Finally, the number of required genes for the calculation of the relevant 
eigengenes is reduced (the tree is 'compacted'). The accuracy of the tree 
is reassessed after removal of each gene. Along the way, several 
self explanatory directories, heatmaps and plots are created and stored under 
\code{saveDir}.
}

\value{
A list with the following components:
\item{call}{The call that created the results.}
\item{wgRes}{A list. The results of WGCNA clustering of the Data by
  \code{\link{wgcna.one.step}}.}
\item{betaRes}{A list. The automatically selected beta (power) parameter
  which was used for the WGCNA clustering. It is the result of the call to
  \code{calculate.beta} using the expression data of \code{mit} conditions(s).}
\item{pigengene}{The pigengene object computed for the clusters, result
  of \code{compute.pigengene}.}
\item{leanrtBn}{A list. The results of \code{\link{learn.bn}} call for 
  learning a Bayesian network using the eigengenes.}
\item{selectedFeatures}{A vector of the names of module eigengenes that
  were considered during the construction of decision trees. 
  If \code{bnNum} >0, this corresponds to the immediate neighbors of the 
  Disease or Effect variable in the consensus network.}
\item{c5treeRes}{A list. The results of \code{\link{make.decision.tree}} call
  for learning decision trees that use the eigengenes as features.}
}

\references{
  Large-scale gene network analysis reveals the significance of extracellular 
  matrix pathway and homeobox genes in acute myeloid leukemia,
  Foroushani A, Agrahari R, Docking R, Karsan A, and Zare H.
  In preparation.}

\author{Amir Foroushani, Habil Zare, Rupesh Agrahari, and  Meghan Short}

\note{
The individual functions are exported to facilitated running the pipeline 
step-by-step in a customized way.}

\seealso{
\code{\link{check.pigengene.input}},
\code{\link{balance}},
\code{\link{calculate.beta}},
\code{\link{wgcna.one.step}},
\code{\link{compute.pigengene}},
\code{\link{learn.bn}}, \code{\link{make.decision.tree}},
\code{\link[WGCNA]{blockwiseModules}}}

\examples{
library(org.Hs.eg.db)
data(aml)
data(mds)
d1 <- rbind(aml,mds)
Labels <- c(rep("AML",nrow(aml)),rep("MDS",nrow(mds)))
names(Labels) <- rownames(d1)
p1 <- one.step.pigengene(Data=d1,saveDir=".", bnNum=10, verbose=1, seed=1, 
      Labels=Labels, toCompact=FALSE, doHeat=FALSE)
plot(p1$c5treeRes$c5Trees[["34"]])
}

\keyword{ optimize }
\keyword{ cluster }
\keyword{ tree }
\keyword{ classif }
\keyword{ models }
