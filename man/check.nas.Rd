\name{check.nas}
\alias{check.nas}
\title{Removes \code{NA}s from a data matrix}

\description{
  Checks \code{Data} for \code{NA} values.}

\usage{
  check.nas(Data, naTolerance=0.05, na.rm=TRUE)
}

\arguments{
  \item{Data}{
    A matrix or data frame containing the expression data, with
    genes corresponding to columns and rows corresponding to
    samples. Rows and columns must be named.}

  \item{naTolerance}{
    A number in the 0-1 range. If the frequency of \code{NA}s in a
    column of \code{Data} is more than this threshold , then that column will be removed.}

  \item{na.rm}{
    If \code{TRUE}, \code{NA}s in the Data will be replaces with the
    average of the column, however, if the frequency of \code{NA}s in
    the column is too high (i.e., more than \code{naTolerance}), the whole column will be removed.}
}

\value{
  A list of:
  \item{cleaned}{The cleaned data with no \code{NA} value. Rows are the
    same as \code{Data}, but some columns may be deleted.}
  \item{tooNaGenes}{A character vector of those genes (i.e., column
    names of \code{Data})  that had too many \code{NA}s, and therefore
    were removed.}
  \item{replacedNaNum}{The number of \code{NA} entries in the matrix that were
    replaced with the average of the corresponding column (gene).}
 }

\author{Habil Zare}

\seealso{
\code{\link{check.pigengene.input}}, \code{\link{Pigengene-package}}
}

\examples{
     data(aml)
     dim(aml)
     aml[1:410]<-NA
     c1 <- check.nas(Data=aml)
     dim(c1$cleaned)
     c1$tooNaGenes
     rm(aml)
}

\keyword{ utilities }
