\name{balance}
\alias{balance}
\title{Balances the number of sampels}

\description{
     Oversamples data by repeating rows such that each condition has
     roughly the same number of samples.}

\usage{
     balance(Data, Labels, amplification = 5, verbose = 0)
}

\arguments{
  \item{Data}{
    A matrix or data frame containing the expression data, with
    genes corresponding to columns and rows corresponding to
    samples. Rows and columns must be named.}
  
  \item{Labels}{
     A (preferably named) vector containing the Labels (condition types) for
     \code{Data}. Names must agree with rows of \code{Data}.}

  \item{amplification}{
    An integer that controls the number of repeats for each condition.
    The number of all samples roughly will be multiplied by this factor
    after oversampling.}

  \item{verbose}{
    The integer level of verbosity. 0 means silent and higher values    
    produce more details of computation.}

}

\value{
  A list of:
  \item{balanced}{The matrix of oversampled data}
  \item{Reptimes}{A vector of integers named by conditions reporting the number of
    repeats for each condition.}
 \item{origSampleInds}{The indices of rows in \code{balanced} that correspond to
    the original samples before oversampling}
}

\author{Habil Zare}

\seealso{
\code{\link{Pigengene-package}},
\code{\link{one.step.pigengene}}, \code{\link{wgcna.one.step}},
\code{\link{compute.pigengene}}}

\examples{
     data(aml)
     data(mds)
     d1 <- rbind(aml,mds)
     Labels <- c(rep("AML",nrow(aml)),rep("MDS",nrow(mds)))
     names(Labels) <- rownames(d1)
     b1 <- balance(Data=d1, Labels=Labels)
     d2 <- b1$balanced
}

\keyword{ utilities }
