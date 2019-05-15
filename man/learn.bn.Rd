\name{learn.bn}
\alias{learn.bn}
\title{Learns a Bayesian network}

\description{
     This function takes as input the eigengenes of all modules and
     learns a Bayesian network using bnlearn package. It builds several
     individual networks from random staring networks by optimizing
     their score. Then, it infers a consensus network from the ones
     with relatively "higher" scores. The default hyper-parameters and
     arguments should be fine for most applications.}

\usage{
learn.bn(pigengene=NULL, Data=NULL, Labels=NULL, bnPath = "bn", bnNum = 100,
  consensusRatio = 1/3, consensusThresh = "Auto", doME0 = FALSE, 
  selectedFeatures = NULL, trainingCases = "All", algo = "hc", scoring = "bde",
  restart = 0, pertFrac = 0.1, doShuffle = TRUE, use.Hartemink = TRUE, 
  bnStartFile = "None", use.Disease = TRUE, use.Effect = FALSE, dummies = NULL,
  tasks = "All", onCluster = !(which.cluster()$cluster == "local"), 
  inds = 1:ceiling(bnNum/perJob), perJob = 2, maxSeconds = 5 * 60, 
  timeJob = "00:10:00", bnCalculationJob = NULL, seed = NULL, verbose = 0,
  naTolerance=0.05)

}
\arguments{
  
  \item{pigengene}{
  An object from \code{\link{pigengene-class}}. 
  The output of \code{\link{compute.pigengene}} function.}

  \item{Data}{
    A matrix or data frame containing the training data with 
    eigengenes corresponding to columns and rows corresponding to 
    samples. Rows and columns must be named.}
  
  \item{Labels}{
    A (preferably named) vector containing the Labels (condition types) for 
    the training data. Names must agree with rows of \code{Data}.}

  \item{bnPath}{
    The path to save the results}
  
  \item{bnNum}{
    The total number of individual networks. In practice, the
    number of learnt networks can be less than \code{bnNum} because 
    some jobs may take too long and be terminated.}

  \item{consensusRatio}{
     A numeric in the range \code{0-1} that determines what
     portion of highly scored networks should be used to build the
     consensus network}

  \item{consensusThresh}{
    A vector of thresholds in the range \code{0-1}. For each 
    threshold \code{t}, a consensus network will be build by
    considering the arcs that are present in at least a fraction
    of \code{t} of the individual networks. Alternatively, if it is
    "Auto" (the default), the threshold will be automatically set
    to the mean plus the standard deviation of the frequencies
    (strengths) of all arcs in the individual networks.}

  \item{doME0}{
    If \code{TRUE}, module 0 (the outliers) will be considered in
    learning the Bayesian network.}

  \item{selectedFeatures}{
    A character vector. If not \code{NULL}, only these features (eigengenes) 
    will be used.}

  \item{trainingCases}{
    A character vector that determines which cases (samples) 
    should be considered for learning the network.}

  \item{algo}{
    The algorithm that bnlean uses for optimizing the score. The
    default is "hc" (hill climbing).  
    See \code{\link[bnlearn]{arc.strength}} for other options and more details.}

  \item{scoring}{
    A character determining the scoring criteria. Use 'bde' and 'bic' for 
    the Bayesian Dirichlet equivalent and Bayesian Information Criterion scores,
    respectively. See \code{\link[bnlearn]{score}} for technical details.}

  \item{restart}{
    The number of random restarts. For technical use only. 
    See \code{\link[bnlearn]{hc}}.}

  \item{pertFrac}{
    A numeric in the range \code{0-1} that determines the number of 
    attempts to randomly insert/remove/reverse an arc on every
    random restart. For technical use only.}

  \item{doShuffle}{
    The ordering of the features (eigengenes) is important in
    making the initial network.  If \code{doShuffle=TRUE}, they will
    be shuffled before making every initial network.}

  \item{use.Hartemink}{
  If \code{TRUE}, Hartemink algorithm will be used to discretize data. Otherwise, interval
  discretization will be applied. See \code{bnlearn:\link[bnlearn:preprocessing]{discretize}.}}

  \item{bnStartFile}{
  Optionally, learning can start from a Bayesian network instead of a random network.
  \code{bnStartFile} should contain a list called \code{selected} and \code{selected$BN} 
  should be an object of \code{\link[bnlearn:bn.class]{bn-class}}.
  Non-technical users can set to \code{"None"} to disable.} 

  \item{use.Disease}{
  If \code{TRUE}, the condition variable \code{Disease} will be included in the network,
  which cannot be the child of any other variable.}

  \item{use.Effect}{
  If \code{TRUE}, the condition variable \code{beAML} will be included in the network,
  which cannot be the parent of any other variable.}

  \item{dummies}{
  A vector of numeric values in the range \code{0-1}. Dummy random variables will be 
  added to  the Bayesian network to check whether the learning process is effective. 
  For development purposes only.}

  \item{tasks}{
  A character vector and a subset of \code{c("learn","harvest","consensus","graph")}
  that identifies the tasks to be done. Useful if part of the analysis was done previously,
  otherwise set to \code{"All"}.}

  \item{onCluster}{
  A Boolean variable that is \code{FALSE} if the learning is not done on a computer cluster.}

  \item{inds}{
  The indices of the jobs that are included in the analysis.}

  \item{perJob}{
  The number of individual networks that are learnt by 1 job.}

  \item{maxSeconds}{
  An integer limiting computation time for each training job that runs locally, 
  i.e., when \code{oncluster=FALSE}.}

  \item{timeJob}{
  The time in \code{"hh:mm:ss"} format requested for each job if they are running on a computer cluster.}

  \item{bnCalculationJob}{
  A script used to submit jobs to the cluster. Set to \code{NULL} if not using a cluster.}

  \item{seed}{
  The random seed that can be set to an integer to reproduce the same results.}

  \item{verbose}{
  Integer level of verbosity. 0 means silent and higher values    
  produce more details of computation.}

  \item{naTolerance}{
   Upper threshold on the fraction of entries per gene that
   can be missing. Genes with a larger fraction of missing
   entries are ignored. For genes with smaller fraction of NA
   entries, the missing values are imputed from their average
   expression in the other samples.  
   See \code{\link{check.pigengene.input}}.}
}

\details{
    For learning a Bayesian network with tens of nodes (eigengenes), \code{bnNum=1000} 
    or higher is recommended. Increasing \code{consensusThresh} generally results 
    in a network with fewer arcs. Nagarajan et al. proposed a fundamental approach that
    determines this hyper-parameter based on the background noise. They use non-parametric 
    bootstrapping, which is not implemented in the current package yet. 

    The default values for the rest of the hyper-parameters should be fine for most applications.}


\value{
A list of: 
  \item{consensusThresh}{The vector of thresholds as described in the arguments.}
  \item{indvPath}{The path where the individual networks were saved.}
  \item{moduleFile}{The file containing data in appropriate format for bnlearn package 
              and the blacklist arcs.}
  \item{scoreFile}{The file containing the record of the successively jobs and the scores of 
                  the corresponding individual networks.}
  \item{consensusFile}{The file containing the consensus network and its BDe and BIC scores.}
  \item{bnModuleRes}{The result of \code{bn.module} function. Useful mostly for development.}
  \item{runs}{A list containing the record of successful jobs.}
  \item{scores}{The list saved in \code{scoreFile}.}
  \item{consensusThreshRes}{The full output of \code{consensus.thresh()} function.}
  \item{consensus1}{The consensus Bayesian network corresponding to the first threshold.
              It is the output of \code{consensus} function and \code{consensus1$BN}
            is an object of \code{\link[bnlearn:bn.class]{bn-class}}.}
  \item{scorePlot}{The output of \code{plot.scores} functions, containing the scores of 
                individual networks.}
  \item{graphs}{The output of \code{plot.graphS} function, containing the BDe score of 
              the consensus network.}
  \item{timeTaken}{An object of \code{\link{difftime}-class} recording the learning wall-time.}
  \item{use.Disease, use.Effect, use.Hartemink}{Some of the input arguments.}
}

\references{
    Hartemink A (2001). Principled Computational Methods for the Validation and
    Discovery of Genetic Regulatory Networks. Ph.D. thesis, School of Electrical 
    Engineering and Computer Science, Massachusetts Institute of Technology.

    Nagarajan, Radhakrishnan, et al. (2010) Functional relationships between genes 
    associated with differentiation potential of aged myogenic progenitors. 
    Frontiers in Physiology 1.}

\author{
Amir Foroushani, Habil Zare, and Rupesh Agrahari}

\note{
Running the jobs on a cluster needs \code{bnCalculationJob} script, which is NOT included in the package yet.}

\seealso{
  \code{\link[bnlearn]{bnlearn-package}}, \code{\link{Pigengene-package}}, 
  \code{\link{compute.pigengene}}
}


\examples{
data(eigengenes33)
ms <- 10:20 ## A subset of modules for quick demonstration
amlE <- eigengenes33$aml[,ms]
mdsE <- eigengenes33$mds[,ms]
eigengenes <- rbind(amlE,mdsE)
Labels <- c(rep("AML",nrow(amlE)),rep("MDS",nrow(mdsE)))
names(Labels) <- rownames(eigengenes)
learnt <- learn.bn(Data=eigengenes, Labels=Labels, 
  bnPath="bnExample", bnNum=10, seed=1)
bn <- learnt$consensus1$BN

## Visualize:
d1 <- draw.bn(BN=bn,nodeFontSize=14)

## What are the children of the Disease node?
childrenD <- bnlearn::children(x=bn, node="Disease")
print(childrenD)

## Fit the parameters of the Bayesian network:
fit <- bnlearn::bn.fit(x=bn, data=learnt$consensus1$Data, method="bayes",iss=10)

## The conditional probability table for a child of the Disease node:
fit[[childrenD[1]]]
}

\keyword{ optimize }
\keyword{ cluster }

