\name{get.used.features}
\alias{get.used.features}

\title{Return the features used in a tree}

\description{
Only some of the features will be automatically selected and used in a 
decision tree. However,  an object of class \code{C5.0} does not have the 
selected feature names explicitly. This function parses the \code{tree}
component and extracts the names of features contributing to the tree.}

\usage{
get.used.features(c5Tree)
}

\arguments{
  \item{c5Tree}{
    A decision tree of class \code{50}}
}

\value{
A character vector of the names of features (module eigengenes)
contributing to the input decision tree.}

\author{Amir Foroushani}

\seealso{
\code{\link{Pigengene-package}},
\code{\link{make.decision.tree}}, \code{\link{compact.tree}}, 
\code{\link{compute.pigengene}}, \code{\link{module.heatmap}},
\code{\link{get.fitted.leaf}}, \code{\link{preds.at}},
\code{\link{Pigengene-package}}
}

\examples{
     ## Data:
     data(aml)
     data(mds)
     data(pigengene)
     d1 <- rbind(aml,mds)

     ## Fiting the trees:
     trees <- make.decision.tree(pigengene=pigengene, Data=d1,
       saveDir="trees", minPerLeaf=15, doHeat=FALSE,verbose=3,
       toCompact=FALSE)
     get.used.features(c5Tree=trees$c5Trees[["15"]])
}

\keyword{ utilities }
\keyword{ tree }
