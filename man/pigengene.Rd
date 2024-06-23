\name{pigengene}
\alias{pigengene}
\docType{data}
\title{An object of class \code{Pigengene}}

\description{
This is a toy example object of class \code{\link{pigengene-class}}.
It is used in examples of \code{\link{Pigengene-package}}.
Gene expression profile of 202  acute myeloid leukemia (AML) cases 
from Mills et al. study. The profile was compared with the profile of 164 
myelodysplastic syndromes (MDS) cases and only the 1000 most differentially expressed
genes are included.}

\usage{data("aml")}

\format{
  An object of \code{\link{pigengene-class}}.}

\value{
  It is an object of \code{\link{pigengene-class}}.}

\details{
    The object is made using \code{\link{compute.pigengene}} function
    from \code{\link{aml}} and \code{\link{mds}} data as shown 
    in the examples. The \code{R CMD build --resave-data} trick was used
    to reduce the size of saved object from 3.1 MB to 1.4 MB.}

\source{\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15061}}

\references{
    Mills, Ken I., et al. (2009). Microarray-based classifiers and prognosis models
    identify subgroups with distinct clinical outcomes and high risk of AML
    transformation of myelodysplastic syndrome. Blood 114.5: 1063-1072.}

\seealso{
  \code{\link{Pigengene-package}}, \code{\link{pigengene-class}},
  \code{\link{one.step.pigengene}},
  \code{\link{mds}}, \code{\link{aml}}, \code{\link{compute.pigengene}},
  \code{\link{project.eigen}}
}

\examples{
library(pheatmap)
data(pigengene)
plot(pigengene,fontsize=12)

## To reproduce:
data(aml)
data(mds)
data(eigengenes33)
d1 <- rbind(aml,mds)
Labels <- c(rep("AML",nrow(aml)),rep("MDS",nrow(mds)))
names(Labels) <- rownames(d1)
modules33 <- eigengenes33$modules[colnames(d1)]
## Computing:
computed <- compute.pigengene(Data=d1, Labels=Labels, modules=modules33,
   saveFile="pigengene.RData", doPlot=FALSE, verbose=3)
class(computed)
plot(computed, fontsize=12, main="Reproduced")
}

\keyword{datasets}
