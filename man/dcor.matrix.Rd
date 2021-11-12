\name{dcor.matrix}
\alias{dcor.matrix}

\title{Computes distance correlation for give matrix}

\description{
This function computes the distance correlation between every pair of columns of
the input data matrix.}

\usage{
dcor.matrix(Data)
}

\arguments{
  \item{Data}{
  A matrix containing the data}
}

\details{
Using for loops, all pairs of columns are passed to \code{link[energy]{dcor}} 
function from \code{link[energy]{energy-package}}.}

\value{
A numeric square matrix. The number of rows and columns is equal to the number
of columns of \code{Data} and they are named accordingly.}

\references{
     Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007), Measuring
     and Testing Dependence by Correlation of Distances, _Annals of
     Statistics_, Vol. 35 No. 6, pp. 2769-2794.

     <URL: http://dx.doi.org/10.1214/009053607000000505>

     Szekely, G.J. and Rizzo, M.L. (2009), Brownian Distance
     Covariance, _Annals of Applied Statistics_, Vol. 3, No. 4,
     1236-1265.

     <URL: http://dx.doi.org/10.1214/09-AOAS312>

     Szekely, G.J. and Rizzo, M.L. (2009), Rejoinder: Brownian Distance
     Covariance, _Annals of Applied Statistics_, Vol. 3, No. 4,
     1303-1308.
}

\author{Habil Zare}

\note{This function uses for loops, which are not efficient for an input matrix
with too many columns.}

\seealso{
\code{link[energy]{dcor}}
}

\examples{
## Data:
data(aml)
dcor1 <- dcor.matrix(Data=aml[,1:5])
dcor1

## Comparison with Pearson:
cor1 <- abs(stats::cor(aml[,1:5]))
## With 202 samples, distance and Pearson correlations do not differ much:
dcor1-cor1 
dcor2 <- dcor.matrix(Data=aml[1:20,1:5])
cor2 <- abs(stats::cor(aml[1:20,1:5]))
## Distance correlation is more robust if fewer samples are available:
dcor2-cor2
plot(dcor2-cor1,cor1-cor2,xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
}

\keyword{ utilities }
