\name{gene.mapping}
\alias{gene.mapping}
\title{Maps gene IDs}

\description{
Takse as input gene IDs in a convention, say REFSEQ, and converts
them to another convention.}

\usage{
gene.mapping(ids, inputType = "REFSEQ", outputType = "SYMBOL",
  leaveNA = TRUE, inputDb = "Human", outputDb = inputDb,
  verbose = 0)
}

\arguments{
  \item{ids}{
  A character vector of input gene IDs}
  
  \item{inputType}{
  The type of input IDs.}

  \item{outputType}{
  The type of output IDs. If it is a character vector, mapping will be
  done for each element.}

  \item{leaveNA}{
  If \code{TRUE},the IDs that were not matched are left with \code{NA}s in
  the second column of the output, otherwise the input IDs are returned.}

  \item{inputDb}{
  The input data base. Use \code{org.Hs.eg.db} for human and
  \code{org.Mm.eg.db} for mouse. The default "Human" character uses the former.}

  \item{outputDb}{
  The output data base. If it is a list, mapping will be done for each element.}

  \item{verbose}{
    The integer level of verbosity. 0 means silent and higher values    
    produce more details of computation.}
}

\details{
It can map homologous genes between species e.g. from mouse to
human.  If more than 1 ID found for an input gene, only one of
them is returned.}

\value{
A matrix of characters with 3 columns: input, output1, and
output2. The last one is guaranteed not to be \code{NA}.}

\references{
Pages H, Carlson M, Falcon S and Li N. AnnotationDbi: Annotation
Database Interface. R package version 1.32.3.
}

\author{
Amir Foroushani, Habil Zare, and Rupesh Agrahari}

\seealso{
\code{\link[AnnotationDbi]{AnnotationDb-class}}, 
\code{\link[org.Hs.eg.db]{org.Hs.eg.db}}
\code{\link[org.Mm.eg.db]{org.Mm.eg.db}}}

\examples{
     library(org.Hs.eg.db)
     g1 <- gene.mapping(ids="NM_001159995")
     print(g1)

     ## Mapping to multiple convention
     library(org.Mm.eg.db)
     g2 <- gene.mapping(ids=c("NM_170730", "NM_001013580"),
        inputType="REFSEQ", inputDb=org.Mm.eg.db,
        outputType=c("SYMBOL","ENTREZID"),
        outputDb=list(org.Hs.eg.db,org.Mm.eg.db), verbose=1)
     print(g2)
}

\keyword{ misc }
