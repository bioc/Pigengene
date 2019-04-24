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
  \code{\link{compute.pigengene}} function. 
  It is a list containing at least the following components:
  \itemize{
    \item{\code{call}}{
      The call that created the results.}
    \item{\code{Reptimes}}{
      A named vector reporting the number of repeats for each condition
      in the oversampling process, which is done by the \code{\link{balance}} function.}
    \item{\code{eigenResults}}{
      The full output of \code{\link[WGCNA]{moduleEigengenes}} function.}
    \item{\code{Data}}{
      The data matrix of gene expression.}
    \item{\code{Labels}}{
      A character vector giving the condition (type) for each sample 
      (row of \code{Data}).} 
    \item{\code{eigengenes}}{
      The matrix of eigengenes ordered based on \code{selectedModules} 
      if provided.}
    \item{\code{membership}}{
      The matrix of weights of genes (rows) in all modules (columns).}
    \item{\code{orderedModules}}{
      The module assignment numeric vector named with genes and ordered
      based on module number.}
    \item{\code{annotation}}{
      A data frame containing labeling information useful in plotting. 
      It has one column named "Condition". Rows have sample names.}
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
