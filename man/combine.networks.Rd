\name{combine.networks}
\alias{combine.networks}

\title{
  Combines two or more networks
}

\description{
  Takes as input two or more adjacency matrices, and the corresponding
  contributions. Computes a combined network (weighted graph) in which the
  weight on an edge between two nodes is an average of the weights on
  the same edge in the input networks.   
}

\usage{
combine.networks(nets, contributions, outPath, midfix="",
    powerVector=1:20, verbose=1, RsquaredCut=0.75, minModuleSize=5,
    doRemoveTOM=TRUE, datExpr, doReturNetworks=FALSE, doSave=FALSE, doIdentifyModule=TRUE)
}

\arguments{
  \item{nets}{
    A list of adjacency matrices (networks), which can be generated using e.g.,
    the \code{WGCNA::\link[WGCNA]{adjacency}} function. Rows and columns must be named.}
  
  \item{contributions}{
    A numeric vector with the same length as nets. In computing the average weight on each
    edge in the combined network, first the edge weights from individual
    networks are multiplied by their corresponding contributions, then the result will be
    divided by the sum of weights of all networks containing this edge.
    }
  
  \item{outPath}{
    A string to the path where plots and results will be saved.}
  
  \item{midfix}{
    An optional string used in the output file names.}
  
  \item{powerVector}{A numeric vector of power values that are tried to
    find the best one.
    See \code{WGCNA::\link[WGCNA]{pickSoftThreshold}} documentation.}

  \item{verbose}{
    The integer level of verbosity. 0 means silent and higher values    
    produce more details of computation.}

  \item{RsquaredCut}{ A threshold in the range [0,1] used to estimate the power. A higher value
  can increase power. For technical use only. See \code{\link[WGCNA]{pickSoftThreshold}} 
  for more details.}

  \item{minModuleSize}{
    The value that controls the minimum number of genes per module.
    See \code{WGCNA::\link[WGCNA]{blockwiseModules}}.}
  
  \item{doRemoveTOM}{
    A boolean determining the big TOM file must remove or not.}  

  \item{datExpr}{
    The expression matrix that
    \code{WGCNA::\link[WGCNA]{blockwiseModules}} uses for fine-tuning
    and removing genes from modules. This is not an ideal behavior by
    WGCNA.}

  \item{doReturNetworks}{
    A boolean value to determine whether to return \code{Network}, which
    is relatively a big matrix (typically GBs). Set to \code{FALSE} not
    to waste memory.
  }
  \item{doSave}{
    A boolean value to determine whether the whole output of this
    function (typically 1-2 GBs) should be saved as \code{combinedNetwork}.
    Set to \code{FALSE} not to waste disk space.
  }
  \item{doIdentifyModule}{
    A boolean value to determine whether modules should be
    identified. Set it to \code{FALSE} if you just need the network, not
    the modules.
  }
}

\value{
  A list with following components
  \item{call}{
    The command that created the results}
  \item{midfix}{
    The input argument}
  \item{Network}{
    The adjacency matrix of the combined network}
  \item{denominators}{
    A matrix, each cell of which is the sum of weights of all networks
    contributing to the edge corresponding to that cell}
  \item{power}{
    The power (beta) value used for the combined network}
  \item{fits}{
    The fit indices calculated for the combined network}
  \item{net}{
    The output of \code{WGCNA::\link[WGCNA]{blockwiseModules}}
    containing the module information in its \code{colors} field}
  \item{modules}{
    The output of \code{WGCNA::\link[WGCNA]{blockwiseModules}}}
  \item{combinedNetworkFile}{
    The path to the saved file containing \code{combinedNetwork}}
}

\note{
    If the networks have different node sets, the combined network will
    be computed on the union of nodes.
}

\seealso{
\code{WGCNA::\link[WGCNA]{blockwiseModules}},
\code{WGCNA::\link[WGCNA]{TOMsimilarity}}, and
\code{WGCNA::\link[WGCNA:pickSoftThreshold]{pickSoftThreshold.fromSimilarity}}
}

\examples{
data(aml)
data(mds)
nets <- list()
## Make the coexpression networks:
nets[["aml"]] <- abs(stats::cor(aml[,1:200]))
nets[["mds"]] <- abs(stats::cor(mds[,1:200]))
## Combine them:
combined <- combine.networks(nets=nets, contributions=c(nrow(aml), nrow(mds)),
                 outPath=".", datExpr=rbind(aml, mds)[,1:200])
print(table(combined$modules))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}% use one of  RShowDoc("KEYWORDS")
\keyword{graphs}% __ONLY ONE__ keyword per line
