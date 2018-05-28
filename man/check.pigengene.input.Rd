\name{check.pigengene.input}
\alias{check.pigengene.input}
\title{Quality check on the imput}

\description{
  Checks \code{Data} and \code{Labels} for \code{NA} values, row and column  
  names, etc.}

\usage{
  check.pigengene.input(Data, Labels, na.rm = FALSE)
}

\arguments{
  \item{Data}{
    A matrix or data frame containing the expression data, with
    genes corresponding to columns and rows corresponding to
    samples. Rows and columns must be named.}
  
  \item{Labels}{
    A (preferably named) vector containing the Labels (condition
    types) for \code{Data}. Names must agree with rows of \code{Data}.}

  \item{na.rm}{
    If \code{TRUE}, \code{NA}s in the Data will be replaces with the
    average of the column, however, if the frequency of \code{NA}s in
    the column is too high, the whole column will be removed.}
}

\value{
   A list of:
  \item{Data}{The checked \code{Data} matrix, \code{NA} possibly removed and rows are
    ordered as names of \code{Labels}.}
  \item{Labels}{The checked vector of \code{Labels}}
}

\author{Habil Zare}

\seealso{
\code{\link{check.nas}}, \code{\link{one.step.pigengene}}, \code{\link{Pigengene-package}}
}

\examples{
     data(aml)
     Labels <- c(rep("AML",nrow(aml)))
     names(Labels) <- rownames(aml)
     c1 <- check.pigengene.input(Data=aml, Labels=Labels,na.rm=TRUE)
     Data <- c1$Data
     Labels <- c1$Labels
}

\keyword{ utilities }
