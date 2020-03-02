\name{message.if}
\alias{message.if}

\title{
  Conditional messaging.}

\description{
  Messages only if verbose is more than 0}

\usage{
message.if(me = NULL, verbose = 0)
}

\arguments{
  \item{me}{
    The Message.}
  
  \item{verbose}{
  A integer.}
}

\value{\code{NULL}}

\author{
  Amir Foroushani}

\examples{
message.if("Hello world!", verbose=1)
}

\keyword{utilities}
