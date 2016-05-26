\name{draw.bn}
\alias{draw.bn}

\title{Draws a Bayesian network}

\description{
Draws the BN using appropriate colors and font size.}

\usage{
draw.bn(BN, plotFile = NULL, inputType = "ENTREZIDat", edgeColor = "blue",
  DiseaseCol = "darkgreen", DiseaseFill = "red", DiseaseChildFill = "pink",
  nodeCol = "darkgreen", nodeFill = "yellow", moduleNamesFile = NULL,
  mainText = NULL, nodeFontSize = 14 * 1.1, verbose = 0)
}

\arguments{
  \item{BN}{
    An object of \code{\link[bnlearn]{bn-class}}}

  \item{plotFile}{
    If not \code{NULL}, the plot will be saved here.}

  \item{inputType}{
    The type of gene IDs in \code{BN}}

  \item{edgeColor}{
    The color of edges}
             
  \item{DiseaseCol}{
    The color of the border of the Disease node}

  \item{DiseaseFill}{
    The color of the area inside the Disease node}

  \item{DiseaseChildFill}{
    The color of the area inside the children of the Disease node}

  \item{nodeCol}{
    The color of the border of the usual nodes excluding Disease 
    and its children}

  \item{nodeFill}{
    The color of the area inside the usual nodes}

  \item{moduleNamesFile}{
     An optional csv file including the information to 
     rename the nodes name.  See code{rename.node}.}

  \item{mainText}{
    The main text shown at the top of the plot}

  \item{nodeFontSize}{
    Adjusts the size of nodes}

  \item{verbose}{
    The integer level of verbosity. 0 means silent and higher values    
    produce more details of computation.} 
}

\value{
A list with following components:
  \item{call}{
    The call that created the results}
  \item{BN}{
    An echo of input \code{BN} argument}
  \item{renamedBN}{
    An object of \code{\link[bnlearn]{bn-class}} when 
    \code{moduleNamesFile} is provided}
  \item{gr}{
    The full output of \code{\link[bnlearn]{graphviz.plot}} function}
}

\author{Habil Zare}

\seealso{
  \code{\link[bnlearn]{bnlearn-package}}, \code{\link{Pigengene-package}}, 
  \code{\link{learn.bn}}, \code{\link[graph]{graph-class}}
}

\examples{
## See lear.bn function.
}

\keyword{ utilities }
