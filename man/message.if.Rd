\name{message.if}
\alias{message.if}

\title{
  Conditional messaging.}

\description{
  Messages only if verbose is more than 0 and write in a text file if provided.}

\usage{
message.if(me=NULL, verbose=0, txtFile=NULL, append=TRUE, ...)
}

\arguments{
  \item{me}{
    The Message. Can be a character vector.}
  
  \item{verbose}{
    A integer.}

  \item{txtFile}{
    The text file in which the message will be written. Set to
    \code{NULL} to disable.}
  
  \item{append}{
    logical. Set to \code{FALSE} to overwrite \code{txtFile}.}
  
  \item{...}{
    Arguments to be passed to \code{capture.output}.}
}

\value{\code{NULL}}

\author{Amir Foroushani}

\examples{
message.if("Hello world!", verbose=1)
}

\keyword{utilities}
