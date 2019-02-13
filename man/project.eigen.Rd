\name{project.eigen}
\alias{project.eigen}
%%
\title{Infers eigengenes for given expression data}

\description{
     This function projects (new) expression data onto the eigengenes
     of modules from another dataset. It is usfull for comparing the
     expression behaviour of modules accross (biologically related yet
     independent) datasets, for evaluating the performance of a
     classifier on new datasets, and for examining the robustness of a
     pattern with regards to missing genes.}

\usage{
project.eigen(Data, saveFile = NULL, pigengene, naTolerance = 0.05, 
  verbose = 0, ignoreModules = c())
}

\arguments{
  \item{Data}{
    A matrix or data frame of expression data to be projected.
    Genes correspond to columns, and rows correspond to samples.
    Rows and columns must be named.  It is OK to miss a few genes
    originally used to compute the eigengenes, thereby,
    projection is robust to choose of platform.}

  \item{saveFile}{
    If not \code{NULL}, where to save the results in \code{.RData} format.}

  \item{pigengene}{
    An object of \code{\link{pigengene-class}}, usually created by  
    \code{\link{compute.pigengene}}}

  \item{naTolerance}{
    Upper threshold on the fraction of entries per gene that
    can be missing. Genes with a larger fraction of missing
    entries are ignored. For genes with smaller fraction of NA
    entries, the missing values are imputed from their average
    expression in the other samples.  
    See \code{\link{check.pigengene.input}}.}

  \item{verbose}{
    The integer level of verbosity. 0 means silent and higher values    
    produce more details of computation.}

  \item{ignoreModules}{
    A vector of integers. In order to speed up the
    projection, it may be desirable to focus only on the
    eigengenes of a few interesting modules. In that case, the
    remaining modules can be listed here and will be ignored
    during projection (Optional).}

}
\details{
    For each module, from the pigengene object, the weight
     (\code{membership}) of each gene is retrieved. The eigengene is
     computed (inferred) on the new data as alinear combination using
     the corresponding weights. The inferred eigengene vector will be
     normalized so that it has the same Euclidean norm as the original
     eigengene vector.}

\value{
A list of:
  \item{projected}{
    The matrix of inferred (projected) eigengenes}
  \item{replacedNaNum}{
    The number of \code{NA} entries in the input \code{Data} that were
    replaced with the the average expression of the corresponding gene}
  \item{tooNaGenes}{
    A character vector of genes that were ignored because they
     had too many \code{NA}s}
  \item{notMatched}{
    A character vector of genes in the original eigengene that
    could not be matched in the given input \code{Data}}
}

\note{
The new data should use the same type of biolocal identifiers
(e.g. Gene Symbols or ENTREZIDs) as the original data for which
the pigengene was constructed. It is, however, not required that
the new data originate from the same type of technology, e.g. the
eigengenes can be based on microarray experiments, whereas the new
data comes from an RNA-Seq experiment. Nor is it necessary that 
the new datset contains measurements for all of the genes from the
original modules.}

\seealso{
\code{\link{Pigengene-package}}, \code{\link{compute.pigengene}}
\code{\link[WGCNA]{moduleEigengenes}}}

\examples{
     ## Data:
     data(aml)
     data(mds)
     data(eigengenes33)
     d1 <- rbind(aml,mds)
     Labels <- c(rep("AML",nrow(aml)),rep("MDS",nrow(mds)))
     names(Labels) <- rownames(d1)
     toyModules <- eigengenes33$modules[colnames(d1)]
     ## Computing:
     p1 <- compute.pigengene(Data=d1, Labels=Labels, modules=toyModules, 
        saveFile="pigengene.RData", doPlot=TRUE, verbose=3)
     ## How robust projecting is?
     p2 <- project.eigen(Data=d1, pigengene = p1, verbose = 1)
     plot(p1$eigengenes[,"ME1"],p2$projected[,"ME1"])
     cor(p1$eigengenes[,"ME1"],p2$projected[,"ME1"])

}

\keyword{ cluster }
\keyword{ classif }
