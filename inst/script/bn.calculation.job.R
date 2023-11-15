## Habil wrote this R script to submit a job to run iterations on cluster.
##-------------------------------------------------------------------------
## Habil added this file from ~/proj/genetwork/code/Rupesh/BN to Pigengene/int/script on 2023-11-14.

library(bnlearn)
library(parallel)
library(getopt)
##library(org.Hs.eg.db)
library(Pigengene)

## Timing:
startTime <- Sys.time()
print("This BN training started at:")
print(startTime)

## Getting the input to the script:
spec = matrix(c(
    'bnInputFile', 'b', 1, "character",
    'moduleNum' , 'm', 1, "character",
    'numClust' , 'n', 1, "integer",
    'perJob' , 'r', 1, "integer",
    'algo', 'a', 1, "character",
    'scoring', 's', 1, "character",
    'indivBNs', 'i', 1, "logical",
    'restart' , 't', 1, "integer",
    'pertFrac' , 'u', 1, "numeric",
    'resultPath', 'e', 1, "character",
    'ind', 'f', 1, "character",
    'bnStartFile', 'l', 1, "character",
    'seed', 'd', 1, "numeric",
    'doShuffle', 'h', 1, "logical",
    'verbose', 'v', 1, "numeric"    
    ), byrow=TRUE, ncol=4)
opt = getopt(spec)
print("The arguments are set as:")
print(opt)


## Loading the inputs:
print(paste("bn.calculate() will load input from: ",opt$bnInputFile,sep=""))
print(paste("And save output in folder:",opt$resultPath,sep=""))
print("The arguments are:")
print(opt)

bns <- bn.calculation(bnInputFile=opt$bnInputFile, moduleNum=opt$moduleNum, numClust=opt$numClust,
                      perJob=opt$perJob, algo=opt$algo, scoring=opt$scoring,
                      indivBNs=opt$indivBNs, restart = opt$restart, pertFrac = opt$pertFrac,
                      resultPath=opt$resultPath, ind=opt$ind, maxSeconds=Inf, seed=opt$seed,
                      verbose=opt$verbose, bnStartFile=opt$bnStartFile, doShuffle=opt$doShuffle)


## Ending:
print("Warnings:")
print(warnings())
print("Time taken:")
print(Sys.time()-startTime)
