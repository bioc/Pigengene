\name{wgcna.one.step}
\alias{wgcna.one.step}

\title{Module identification}

\description{
     This function is a wrapper function for
     WGCNA::\code{\link[WGCNA]{blockwiseModules}} and passes its arguments to it. 
     Some other arguments are fixed.}

\usage{
wgcna.one.step(Data, power, saveDir=".", blockSize = "All", saveTOMs = FALSE, 
   doThreads=FALSE, verbose = 0, seed = NULL)
}

\arguments{

  \item{Data}{
    A matrix or data frame containing the expression data, with
    genes corresponding to columns and rows corresponding to
    samples. Rows and columns must be named.}

  \item{power}{
    Soft-thresholding power for network construction}

   \item{saveDir}{
     The directory to save the results and plots. \code{NULL} will
     disable saving.}

   \item{blockSize}{
     The size of block when the data is too big. If not "All"
     (default) may introduce artifacts.}

   \item{saveTOMs}{
     Boolean determining if the TOM data should be saved, which
     can be hundreds of MBs and useful for identifying hubs.}

  \item{doThreads}{
    Boolean. Allows WGCNA to run a little faster using multi-threading
    but might not work on all systems.}

  \item{verbose}{
    The integer level of verbosity. 0 means silent and higher values    
    produce more details of computation.}

  \item{seed}{
    Random seed to ensure reproducibility.}
}

\details{
\code{Data, power, blockSize, saveTOMs, verbose,} and \code{seed}
are passed to WGCNA::\code{\link[WGCNA]{blockwiseModules}}.
}


\value{
A list with following components
\item{call}{
    The command that created the results}
\item{genes}{
    The names of \code{Data} columns}
\item{modules}{
    A numeric vector, named by \code{genes}, that reports the module
  (clustering) assignments.}

\item{moduleColors}{
    A character vector, named by \code{genes}, that reports the
  color of each gene according to its module assignment}
\item{net}{
    The full output of \code{\link[WGCNA]{blockwiseModules}} function}
\item{netFile}{
    The file in which the net object is saved}
\item{power}{
    An echo of the \code{power} argument.}
}

\references{
Langfelder P and Horvath S, WGCNA: an R package for weighted
correlation network analysis. BMC Bioinformatics 2008, 9:559}

\seealso{
  \code{\link[WGCNA]{blockwiseModules}},
  \code{\link[WGCNA]{pickSoftThreshold}},
  \code{\link{calculate.beta}}
}

\examples{
     data(aml)
     wgRes <- wgcna.one.step(Data=aml[,1:200], seed=1, power=7,
                                 saveDir="wgcna", verbose=1)
}

\keyword{ cluster }
