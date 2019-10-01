\name{save.if}
\alias{save.if}

\title{
  Saves an object verbosely.
}

\description{
  Saves an R object, and reports the size of the saved object in
  memory and on file.
}

\usage{
save.if(x1, file, verbose = 1)
}

\arguments{
  \item{x1}{
    The object to be saved.}
  
  \item{file}{
    Where to save. If \code{NULL}, nothing will be saved.}
  
  \item{verbose}{
  A numeric determining how much detail will be printed.}
}
\value{
  A list including file, and a vector of sizes of the object in memory and on file.}

\author{
  Amir Foroushani, and Habil Zare}

\seealso{
\code{\link{message.if}}
}

\examples{
m1 <- matrix(0, nrow=1000, ncol=1000)
save.if(m1, file="./m1.RData", verbose=3)
}

\keyword{utilities}
