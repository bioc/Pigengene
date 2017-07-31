\name{pheatmap.type}
\alias{pheatmap.type}

\title{
Plots heatmap with clustering only within types.}

\description{
     This function first performs hierarchical clustering on samples
     (rows of data) within each condition. Then, plots a heatmap
     without further clustering of rows.}

\usage{
pheatmap.type(Data, annRow, type = colnames(annRow)[1], doTranspose=FALSE, ...)
}

\arguments{
  \item{Data}{
    A matrix with samples on rows and features (genes) on columns.}

  \item{annRow}{
    A data frame with 1 column or more. Row names must be the same 
    as row names of Data.}

  \item{type}{
    The column of \code{annRow} used for determining the condition}

  \item{doTranspose}{
    If \code{TRUE}, the matrix will be transposed for the final
    plot. E.g., if the genes are on the columns of \code{Data}, they
    will be shown on rows of the heatmap.}

  \item{\dots}{
    Additional arguments passed to \code{pheatmap} function.}
}

\value{
A list of:
  \item{pheatmapS}{
    The results of pheatmap function for each condition}
  \item{pheat}{
    The output of final pheatmap function applied on all data}
  \item{ordering}{
    The ordering of the rows in the final heatmap}
  \item{annRowAll}{
    The row annotation used in the final heatmap}
}

\author{Habil Zare}

\note{
 If \code{type} is not determined, by default the first column of
 \code{annRow} is used.}

\seealso{
\code{\link{eigengenes33}}}

\examples{
     data(eigengenes33)
     d1 <- eigengenes33$aml
     d2 <- eigengenes33$mds
     Disease <- c(rep("AML",nrow(d1)), rep("MDS",nrow(d2)))
     Disease <- as.data.frame(Disease)
     rownames(Disease) <- c(rownames(d1), rownames(d2))
     p1 <- pheatmap.type(rbind(d1,d2),annRow=Disease,show_rownames=FALSE)
}

\keyword{ hplot }
\keyword{ cluster }
