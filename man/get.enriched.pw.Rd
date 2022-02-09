\name{get.enriched.pw}
\alias{get.enriched.pw}
\title{Performs pathway over representation analysis}

\description{
Takes as input a vector or list of gene IDs in any convention, and performs
over representation analysis.}

\usage{
get.enriched.pw(genes, idType, pathwayDb, ont=c("BP", "MF", "CC"),
                Org="Human", OrgDb=NULL, outPath, pvalueCutoff=0.05, 
                pAdjustMethod="BH", fontSize=14, verbose=0)
}

\arguments{
  \item{genes}{
  A character vector or a named list of genes for which pathway over representation 
  analysis to be done.}
  
  \item{idType}{
  A string describing the type of input gene ID e.g., "ENTREZID", "REFSEQ", 
  "SYMBOL".}

  \item{pathwayDb}{
  A character vector determining which enrichment database to be used e.g., "GO", 
  "KEGG", "REACTOME", or "NCG".}

  \item{ont}{
  GO ontology terms to be analysed e.g., "BP", "MF" or "CC". Default is all three.}

  \item{Org}{
  A character string equal to "Human" or "Mouse" determining the reference
  organism to be used. For "Human" and "Mouse" \code{org.Hg.eg.db} and 
  \code{org.Mm.eg.db} will be used, respectively. If \code{Org} is not
  \code{NULL}, \code{OrgDb} must be \code{NULL}.}

  \item{OrgDb}{
  The reference data base to be used. Use e.g. \code{org.Ce.eg.db} for 'Celegans'
  when analysing Celegans data. If \code{OrgDb} is not \code{NULL}, \code{Org} 
  must be \code{NULL}.}

  \item{outPath}{
  A file path where results will be saved.}

  \item{pvalueCutoff}{
  A numerical value that determines a cutoff of adjusted pValue.}

  \item{pAdjustMethod}{
  A string passed to the clusterProfiler::\code{\link[clusterProfiler]{enrichGO}} 
  function to determine the method for adjusting the p-value. Options include
  "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".}

  \item{fontSize}{
  A numerical value that determines the font size of the y-axis and the title
  in the plot.}


  \item{verbose}{
  The integer level of verbosity. 0 means silent and higher values    
  produce more details of computation.}
}

%% \details{}

\value{
A list:
  \item{EnrichResults}{A list of output of enrichment analysis for different
                       database analyzed.}
  \item{noEnrichment}{A vector of database names in which no enriched pathways 
                      were found.}
The output is saved for each selected module under the ``moduleName\_enrichment'' 
folder. There is a subfolder that includes an excel file and plot(s). Each sheet 
in the excel file corresponds to a pathway database (KEGG in the below example). 
Each row is an overrepresented pathway.
}

\references{
Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He,
clusterProfiler: an R package for comparing biological themes among
gene clusters. OMICS: A Journal of Integrative Biology 2012, 16(5):284-287

Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for
reactome pathway analysis and visualization. Molecular BioSystems 2016, 12(2):477-479
}

\author{
Isha Mehta, and Habil Zare}

\seealso{
\code{\link[clusterProfiler]{enrichGO}},
\code{\link[clusterProfiler]{enrichKEGG}},
\code{\link[DOSE]{enrichNCG}}, \code{\link[ReactomePA]{enrichPathway}}
}

\examples{
library(org.Hs.eg.db)
genes <- c("NM_170730", "NM_001013580", "NM_002142", "NM_003417", "NM_000082",
           "NM_006158", "NM_006047", "NM_022356", "NM_003979", "NM_001030", "NM_022872")
p1 <- get.enriched.pw(genes=genes, idType="REFSEQ", pathwayDb="GO", Org="Human", 
                      outPath=getwd(), verbose=1)
}

\keyword{ misc }
