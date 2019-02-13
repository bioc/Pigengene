\name{calculate.beta}
\alias{calculate.beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates an appropriate power value}

\description{
The \code{\link[WGCNA:pickSoftThreshold]{WGCNA}} package assumes that in the coexpression network 
the genes are connected with a power-law distribution. Therefore, it need a 
soft-thresholding power for network construction, which is estimated 
by this auxiliary function.}

\usage{
calculate.beta(saveFile = NULL, RsquaredCut = 0.8, Data,  doThreads=FALSE, 
  verbose = 0)
}

\arguments{
  \item{saveFile}{
  The file to save the results in. Set to \code{NULL} to disable.}

  \item{RsquaredCut}{
  A threshold in the range [0,1] used to estimate the power. A higher value
  can increase power. For technical use only. See \code{\link[WGCNA]{pickSoftThreshold}} 
  for more details.}

  \item{Data}{
  A matrix or data frame containing the expression data, with genes 
  corresponding to columns and rows corresponding to samples. Rows and columns 
  must be named.}

  \item{doThreads}{
    Boolean. Allows WGCNA to run a little faster using multi-threading
    but might not work on all systems.}

  \item{verbose}{
  The integer level of verbosity. 0 means silent and higher values    
  produce more details of computation.} 

}

\value{
A list of:
    \item{sft}{The full output of \code{\link[WGCNA]{pickSoftThreshold}}
    function}
    \item{power}{The estimated power (beta) value}
    \item{powers}{The numeric vector of all tried powers}
    \item{RsquaredCut}{The value of input argument \code{RsquaredCut}}
}

\references{
Langfelder P and Horvath S, WGCNA: an R package for weighted correlation network
analysis. BMC Bioinformatics 2008, 9:559}

\seealso{
\code{\link[WGCNA]{pickSoftThreshold}},
\code{\link[WGCNA]{blockwiseModules}},
\code{\link{one.step.pigengene}},
\code{\link{wgcna.one.step}}}

\examples{
     data(aml)
     p1 <- calculate.beta(Data=aml[,1:200])
}

\keyword{ cluster }

