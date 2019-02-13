\name{plot.pigengene}
\alias{plot.pigengene}

\title{Plots and saves a \code{pigengene} object}

\description{
Plots a couple of heatmaps of expression of the eigengenes,
weights (memberships), and so on. Saves the plots in png format.}


\usage{
\method{plot}{pigengene}(x, saveDir = NULL, 
  DiseaseColors = c("red", "cyan"), 
  fontsize = 35, doShowColnames = TRUE, fontsizeCol = 25, 
  doClusterCols = ncol(pigengene$eigengenes) > 1, 
  verbose = 2, doShowRownames = "Auto", 
  pngfactor = max(2, ncol(pigengene$eigengenes)/16), ...)
}

\arguments{
  \item{x}{
    The object from \code{\link{pigengene-class}} computed by 
    \code{\link{compute.pigengene}}.}
  
  \item{saveDir}{
    The dirctory for saving the plots}

  \item{DiseaseColors}{
    A vector of characters determining color for each disease}

  \item{fontsize}{
  Passd to \code{\link{pheatmap.type}}}

  \item{doShowColnames}{
  Boolean}

  \item{fontsizeCol}{
  Numeric}

  \item{doClusterCols}{
  Boolean}

  \item{verbose}{
    The integer level of verbosity. 0 means silent and higher values    
    produce more details of computation.}

  \item{doShowRownames}{
    Boolean}
  
  \item{pngfactor}{
    A numeric adjusting the size of the png files}

  \item{\dots}{
    Passd to \code{pheatmap.type} function}

}

\details{
Many of the arguments are passed to \code{\link[pheatmap]{pheatmap}}.}

\value{
A list of:
  \item{heat}{The full output of \code{pheatmap} functionion}
  \item{heatNotRows}{The full output of \code{pheatmap.type} function}
}

\references{
     Large-scale gene network analysis reveals the significance of
     extracellular matrix pathway and homeobox genes in acute myeloid
     leukemia, Foroushani A, Agrahari R, Docking R, Karsan A, and Zare
     H.  In preparation.}

\author{Habil Zare ad Amir Foroushani}

\seealso{
\code{\link{Pigengene-package}},
\code{\link{compute.pigengene}},
\code{\link{pheatmap.type}}}

\examples{
     ## Data:
     data(aml)
     data(mds)
     data(eigengenes33)
     d1 <- rbind(aml,mds)
     Labels <- c(rep("AML",nrow(aml)),rep("MDS",nrow(mds)))
     names(Labels) <- rownames(d1)
     Labels <- c(rep("AML",nrow(eigengenes33$aml)),rep("MDS",nrow(eigengenes33$mds)))
     names(Labels) <- rownames(d1)
     toyModules <- eigengenes33$modules[colnames(d1)]
     ## Computing:
     p1 <- compute.pigengene(Data=d1, Labels=Labels, modules=toyModules, 
        saveFile="pigengene.RData", doPlot=TRUE, verbose=3)
     plot(p1,saveDir="plots")

}


\keyword{ cluster }
