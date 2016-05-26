\name{eigengenes33}
\alias{eigengenes33}
\docType{data}
\title{
Eigengenes of 33 modules}
\description{
This list contains partial eigengenes computed from  
AML and MDS gene expression profiles provided by Mills et al. 
These data are included to illustrate how to use \code{\link{Pigengene-package}}
and also to facilitate reproducing the results presented in the corresponding 
paper.}

\usage{data(eigengenes33)}
\format{
  A list}

\value{
  It is a list of 3 objects:

 \code{aml}
 {A 202 by 34 matrix. 
 Each column reports the values of a module eigengene for AML cases.}

 \code{mds}
 {A 164 by 34 matrix for MDS cases with columns similar to aml.}
 
 \code{modules} 
 {A numeric vector of length 9166 labeling members of each module. Named by Entrez ID.}
}

\details{
    The top 9166 differentially expressed genes were identified and their expressions
    in AML were used for identifying 33 modules. The first column, ME0, corresponds 
    to module 0 (outliers) and is usually ignored. The eigengene for each module was
    obtained using \code{\link{compute.pigengene}} function. Oversampling 
    was performed with \code{amplification=5} to adjust for unbalanced sample-size.}

\source{\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15061}}

\references{
    Mills, Ken I., et al. (2009). Microarray-based classifiers and prognosis models
    identify subgroups with distinct clinical outcomes and high risk of AML
    transformation of myelodysplastic syndrome. Blood 114.5: 1063-1072.}

\seealso{
  \code{\link{Pigengene-package}}, \code{\link{compute.pigengene}}, 
  \code{\link{aml}}, \code{\link{mds}},
  \code{\link{learn.bn}}
}

\examples{
library(pheatmap)
data(eigengenes33)
pheatmap(eigengenes33$aml,show_rownames=FALSE)
## See Pigengene::learn.bn() documentation for more examples.
}

\keyword{datasets}
