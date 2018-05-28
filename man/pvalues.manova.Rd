\name{pvalues.manova}
\alias{pvalues.manova}

\title{Computes pvalues for multi-class differential expression}

\description{
Passes the arguments to \code{\link[stats]{manova}}, which performs multi-class
analysis of variance.}

\usage{
pvalues.manova(Data, Labels)
}

\arguments{
  \item{Data}{
    A matrix or data frame containing the (expression) data, with
    genes corresponding to columns and rows corresponding to
    samples. Rows and columns must be named.}
  
  \item{Labels}{
    A (preferably named) vector containing the Labels (condition
    types).  Names must agree with rows of \code{Data}}
}

\value{
A list with following elements:
  \item{call}{
    The call that created the results}
  \item{pvals}{
    The matrix of pvalues with columns "pValue", "FDR", "Bonferroni".  
    Rows are named according to genes, the columns of \code{Data}.}
  \item{manovaFit}{
    The full output of \code{\link[stats]{manova}} function.}
}

\references{
     Krzanowski, W. J. (1988) _Principles of Multivariate Analysis. A
     User's Perspective._ Oxford.

     Hand, D. J. and Taylor, C. C.  (1987) _Multivariate Analysis of
     Variance and Repeated Measures._ Chapman and Hall.

     B. L. Welch (1951), On the comparison of several mean values: an
     alternative approach.
}

\author{Amir Foroushani}

\note{
\code{\link[stats]{oneway.test}} function is a better generalizatoion to 
Welch's t-tst from 2-calsses to multi-class because it dose not assume that 
the variaces are necessarly equal. However, in practice, with "enough number 
of samples", the two approaches will lead to similar p-values.}

\seealso{
\code{\link[stats]{oneway.test}}, \code{\link[stats]{manova}}, 
\code{\link{compute.pigengene}}}

\examples{
     data(eigengenes33)
     d1 <- rbind(eigengenes33$aml,eigengenes33$mds)
     Labels <- c(rep("AML",nrow(eigengenes33$aml)),rep("MDS",nrow(eigengenes33$mds)))
     names(Labels) <- rownames(d1)
     ps <- pvalues.manova(Data=d1, Labels=Labels)
     plot(log10(ps$pvals[,"Bonferroni"]))
}

\keyword{ utilities }
