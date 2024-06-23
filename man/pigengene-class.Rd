\name{pigengene-class}
\docType{class}

% Classes
\alias{pigengene-class}

% Standard generic methods:
\alias{plot, pigengene-method}


\title{The pigengene class}

\description{
A pigengene object holds the eigengenes, weights (memberships) 
and other related information.}

\details{
  A object of class \code{pigengene} is the output of 
  \code{\link{compute.pigengene}} function. It is a list.
  If \code{doMinimize=TRUE}, only the minimal elements needed to project
  eigengenes in a new dataset are included
  (i.e., see \code{project.eigen(pigengene=NULL)[["projectionaries"]]}).
  Otherwise, it contains at least the following components:
  \itemize{
    \item{\code{call}}{
      The call that created the results.}
    \item{\code{Reptimes}}{
      A named vector reporting the number of repeats for each condition
      in the oversampling process, which is done by the \code{\link{balance}} function.}
    \item{\code{eigenResults}}{
      A list including at least \code{eigengenes} and
      \code{varExplained}. If \code{doWgcna=TRUE}, then 
      this list will be the full output of
      \code{\link[WGCNA]{moduleEigengenes}} function with some fixes, e.g.,
      we change \code{eigengenes} to a matrix,
      and use genes as its row names. Also, \code{varExplained} is named
      according to modules. Setting \code{doWgcna=TRUE} leads to more
      memory usage and a larger \code{Pigengene} object likely, with no
      advantage.}
    \item{\code{Data}}{
      The data matrix of gene expression.}
    \item{\code{Labels}}{
      A character vector giving the condition (type) for each sample 
      (row of \code{Data}).} 
    \item{\code{eigengenes}}{
      The matrix of eigengenes ordered based on \code{selectedModules} 
      if provided. Rows correspond to samples.}
    \item{\code{membership}}{
      The matrix of weights of genes (rows) in all modules (columns).}
    \item{\code{orderedModules}}{
      The module assignment numeric vector named with genes and ordered
      based on module number.}
    \item{\code{annotation}}{
      A data frame containing labeling information useful in plotting. 
      It has a column named "Condition". Rows have sample names.}
    \item{\code{saveFile}}{
      The file where the \code{pigengene} object is saved.}
    \item{\code{weightsCsvFile}}{
      The file containing the weights in csv format. See \code{dOrderByW=TRUE}.}
     \item{\code{weights}}{
      The weight matrix, which is also saved in csv format. It has
      more columns than \code{membership} but rows may be in a different
      order if \code{dOrderByW=TRUE}.}
     \item{\code{heavyToLow}}{
       If \code{dOrderByW=TRUE}, this will be the ordering of genes
       according to the modules the belong to, where the genes in each
       module are ordered based on the absolute value of the weights in
       that module. Also, the genes in the csv file are in this order.}
   }

For 2 or more groups (conditions), additional (optional) components include:
  \itemize{
    \item{\code{pvalues}}{
      A numeric matrix with columns "pValue", "FDR", and "Bonferroni".
      Rows correspond to modules. The null hypothesis is that the eigengene is
      expressed with the same distribution in all groups (conditions).}
    \item{\code{log.pvalues}}{
      A data frame with 1 column containing the logarithm of 
      Bonferroni-adjusted pvalues in base 10.}
    }
}


\seealso{
\code{\link{Pigengene-package}},  \code{\link{plot.pigengene}},
\code{\link{wgcna.one.step}}, \code{\link{compute.pigengene}}, 
\code{\link{learn.bn}}, \code{\link{make.decision.tree}}}

\keyword{methods}
\keyword{classes}
