\name{mds}
\alias{mds}
\docType{data}
\title{MDS gene expression profile}

\description{
Gene expression profile of 164 myelodysplastic syndromes (MDS) cases 
from Mills et al. study. The profile was compared with the profile of 
202  acute myeloid leukemia (AML) cases and only the 1000 most differentially expressed 
genes are included.}

\usage{data("mds")}

\format{
  A numeric matrix}

\value{
  It is a \code{164*1000} numeric matrix.}

\details{
    The columns and rows are named according to the genes Entrez, 
    and patient IDs, respectively. The original data was produced using 
    Affymetrix Human Genome U133 Plus 2.0 Miccoaray.Mills et al. study
    is part of the MILE Study (Microarray Innovations In LEukemia) program,
    and aimed at prediction of AML transformation in MDS.}

\source{\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15061}}

\references{
    Mills, Ken I., et al. (2009). Microarray-based classifiers and prognosis models
    identify subgroups with distinct clinical outcomes and high risk of AML
    transformation of myelodysplastic syndrome. Blood 114.5: 1063-1072.}
\note{
This profile includes data of the 25 chronic myelomonocytic leukemia (CMLL) cases
that can have different expression signatures according to Mills et al.}

\seealso{
  \code{\link[Pigengene]{Pigengene-package}}, 
  \code{\link[Pigengene]{one.step.pigengene}},
  \code{\link[Pigengene]{aml}}, \code{\link{compute.pigengene}}
}

\examples{
library(pheatmap)
data(mds)
pheatmap(mds[,1:20],show_rownames=FALSE)
}

\keyword{datasets}
