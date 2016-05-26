\name{aml}
\alias{aml}
\docType{data}
\title{AML gene expression profile}

\description{
Gene expression profile of 202  acute myeloid leukemia (AML) cases 
from Mills et al. study. The profile was compared with the profile of 164 
myelodysplastic syndromes (MDS) cases and only the 1000 most differentially expressed
genes are included.}

\usage{data("aml")}

\format{
  A numeric matrix}

\value{
  It is a \code{202*1000} numeric matrix.}

\details{
    The columns and rows are named according to the genes Entrez, 
    and patient IDs, respectively. The original data was produced using 
    Affymetrix Human Genome U133 Plus 2.0 Miccoaray. Mills et al. study
    is part of the MILE Study (Microarray Innovations In LEukemia) program,
    and aimed at prediction of AML transformation in MDS.}

\source{\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15061}}


\references{
    Mills, Ken I., et al. (2009). Microarray-based classifiers and prognosis models
    identify subgroups with distinct clinical outcomes and high risk of AML
    transformation of myelodysplastic syndrome. Blood 114.5: 1063-1072.}

\seealso{
  \code{\link{Pigengene-package}}, 
  \code{\link{one.step.pigengene}},
  \code{\link{mds}}, \code{\link{pigengene}}
}

\examples{
library(pheatmap)
data(aml)
pheatmap(aml[,1:20],show_rownames=FALSE)
}

\keyword{datasets}
