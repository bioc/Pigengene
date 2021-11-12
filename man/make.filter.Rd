\name{make.filter}
\alias{make.filter}

\title{
  Computes the filter based on a similarity network}

\description{
    Takes as input the similarity matrix of a graph (i.e., \code{network}) and an \code{epsilon} value.
    It computes a filter graph using the \code{epsilon} threshold.
    The dimention of the output filter matrix is the same as the input similarity network.
    It also produces two plots showing the weighted degrees of the input graph and degrees of the filter,
    respectively.
}

\usage{
make.filter(network, epsilon, outPath=NULL)
}

\arguments{
  \item{network}{
    A matrix of similarity for the network.}
  
  \item{epsilon}{
    A threshold for deciding which edges to keep. If the similarity is less than 1/\code{epsilon}
    (i.e., distance > \code{epsilon}), the edge will be removed, and it will be kept in the filter graph
    otherwise. }

  \item{outPath}{
    A string determining the path where plots and results will be saved.}
}


\value{
A list with the following components:
\item{filt}{A matrix representing adjacency matrix of the computed filter graph. If the distance
between two nodes in the similarity matrix is higher than epsilon, those nodes are connected in
the filter graph (i.e., the corresponding entry in the adjacency matrix is 1).
Otherwise, the corresponding entry is 0.}

\item{epsilon}{The \code{epsilon} input.}
}

\author{Habil Zare and Neda Emami.}

\seealso{
\code{\link{one.step.pigengene}},
\code{\link{apply.filter}}
}

\examples{
data(aml)
data(mds)
d1 <- rbind(aml,mds)[, 1:200]
Labels <- c(rep("AML",nrow(aml)),rep("MDS",nrow(mds)))
names(Labels) <- rownames(d1)

p0 <- one.step.pigengene(Data=d1, saveDir=".", verbose=1,
                         seed=1, Labels=Labels, naTolerance=0.5,
			 RsquaredCut=0.8, doNetOnly=TRUE)

##making the filter
made <- make.filter(network=p0$Network, epsilon=0.7, outPath=".")
}

\keyword{ optimize }
