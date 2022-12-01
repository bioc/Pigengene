\name{determine.modules}
\alias{determine.modules}

\title{
  Identifies modules of the network
}

\description{
    Takes as input a \code{network} (i.e., weighted graph) and identifies modules
    (i.e., clusters of similar genes) using  \code{WGCNA::\link[WGCNA]{blockwiseModules}}.
    It also produces a plot showing the number of genes in each module.
}

\usage{
determine.modules(network, outPath, midfix="", powerVector=1:20,
		verbose=1, RsquaredCut=0.75, minModuleSize=5, doRemoveTOM=FALSE,
                datExpr, doSave=FALSE)
}

\arguments{
  \item{network}{
    An adjacency matrix of the network that is built using \code{\link{combine.networks}}.}

  \item{outPath}{
    A string to the path where plots and results will be saved.}
  
  \item{midfix}{
    An optional string used in the output file names.}
  
  \item{powerVector}{
    A numeric vector of integer values that are tried to find the best power.
    See \code{WGCNA::\link[WGCNA]{pickSoftThreshold}}.}

  \item{verbose}{
    The integer level of verbosity, where 0 means silent and higher values    
    produce more details.}

  \item{RsquaredCut}{
    A threshold in the range [0,1] used to estimate the power. A higher value
    can increase power. For technical use only. See \code{\link[WGCNA]{pickSoftThreshold}} 
    for more details.}

  \item{minModuleSize}{
    The value that controls the minimum number of genes per module.
    See \code{WGCNA::\link[WGCNA]{blockwiseModules}}.}
  
  \item{doRemoveTOM}{
    A boolean determining whether the big TOM file must remove or not.}

  \item{datExpr}{
    The expression matrix that \code{WGCNA::\link[WGCNA]{blockwiseModules}}
    uses for fine-tuning and removing genes from modules. This is not an ideal behavior by
    WGCNA.}

  \item{doSave}{
    A boolean value to determine whether the whole output of this
    function (typically 1-2 GBs) should be saved as \code{combinedNetwork}.
    Set to \code{FALSE} not to waste disk space.}
}

\value{
A list with the following components:
\item{call}{The call that created the results.}
\item{midfix}{The \code{midfix} input.}
\item{power}{The integer value of the estimated power computed by
\code{\link[WGCNA]{pickSoftThreshold.fromSimilarity}}.}
\item{fits}{The \code{fitIndices} output from \code{\link[WGCNA]{pickSoftThreshold.fromSimilarity}}.}
\item{modules}{A vector that representing the identified modules. Its length is equal to the number
of nodes in the network, named by node names (i.e., row names of \code{network}), and values are the
corresponding module numbers.}
\item{net}{The full output of the \code{\link[WGCNA]{blockwiseModules}} function.}
}

\author{Neda Emami and Habil Zare.}


\seealso{
\code{\link{apply.filter}},
\code{\link{combine.networks}},
\code{\link{make.filter}}
}

\examples{
data(aml)

##Making the coexpression network
network <- abs(stats::cor(aml[,1:200]))

##Identifying modules
identifiedMod <- determine.modules(network=network, outPath=".", datExpr=aml[,1:200])
print(table(identifiedMod$modules))
}

\keyword{ graph }
\keyword{ cluster }
