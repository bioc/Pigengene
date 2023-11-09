\name{apply.filter}
\alias{apply.filter}

\title{
  Applies a given filter on the data
}

\description{
    Takes as input \code{gamma} and \code{epsilon} values and a filter graph, which is represented by
    an adjacency matrix named \code{filt}. Applies the filter on the data in either of the two ways:
    a) with normalization of the filter by degrees in the graph,
    b) without normalization.
}

\usage{
apply.filter(gamma, filt, Data, doNormalize=FALSE)
}

\arguments{
  \item{gamma}{
    This value is in the [0,1] range and determines the weight of the filter data.
    Setting to 0 will result in not filtering at all.}
  
  \item{filt}{
  It is a binary matrix computed by the \code{\link{make.filter}} function.}

  \item{Data}{
  A matrix or data frame (or list of matrices or data frames) containing
  the expression data, with genes corresponding to columns and
  rows corresponding to samples. Rows and columns 
  must be named. For example, for RNA-Seq data, log(RPKM+1) can be used.}

  \item{doNormalize}{
  If \code{TRUE}, the filter will be normalized by the degree in the graph using
  the filt * D^(-1), where D is a diagonal matrix with degrees of filt on its diagonal.}
}


\value{
\item{filtered}{A filtered matrix computed using the gamma*sData %*% filtN + (1-gamma)*sData
formula, where sData is the scaled \code{Data} and filtN is the normalized or unormalized filter.}
}


\author{Habil Zare and Neda Emami.}

\seealso{
\code{\link{make.filter}},
\code{\link{determine.modules}}
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

##Making the filter
made <- make.filter(network=p0$Network, epsilon=0.7, outPath=".")

##Applying the filter
f1 <- apply.filter(gamma=0.5, filt=made$filt, Data=d1)
}

\keyword{ optimize }

