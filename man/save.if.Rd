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
save.if(x1, file, compress=TRUE, verbose=1, ...)
}

\arguments{
  \item{x1}{
    The object to be saved.}
  
  \item{file}{
    Where to save. If \code{NULL}, nothing will be saved.}

  \item{compress}{
    A Boolean or character sent to the \code{\link{save}} function. The default
    \code{TRUE} leads to compression using gzip. With "xz", maximum
    compression is obtained in expense of more save and load time.}

  \item{verbose}{
    A numeric determining how much detail will be printed.}

  \item{...}{
    Optional arguments to be passed to the \code{\link{save}} function.}
}
\value{
  A list including file, and a vector of sizes of the object in memory and on file.}

\author{
  Amir Foroushani, and Habil Zare}

\seealso{
\code{\link{message.if}}, \code{\link{save}}
}

\examples{
m1 <- matrix(0, nrow=1000, ncol=1000)
save.if(m1, file="./m1.RData", verbose=3)
}

\keyword{utilities}
