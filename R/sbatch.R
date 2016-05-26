sbatch <- function(
    tempFolder="~/temp/", sourceFile, namePrefix=basename(sourceFile), 
    wait=1, arguments=c(), sourceBashrcBy="~/.bashrc", 
    doSubmit=TRUE, doDeleteTempFile=FALSE, jobName=NULL, 
    ste=NULL, sto=NULL, doTalk=TRUE, queue="normal", 
    time1="47:00:00", maxJobs) 
{
    startTime <- Sys.time()
    tempFolder <- paste(tempFolder, "/sbatch/", basename(sourceFile), 
                        "/", sep="")
    dir.create(tempFolder, recursive=TRUE, showWarnings=FALSE)
    argumentsString <- paste(arguments, collapse="_")
    if (is.null(jobName)) 
        jobName <- paste(namePrefix, argumentsString, sep="_")
    tempFile <- paste(tempFolder, "/", jobName, "_", startTime, 
                      "_sbatch.mpi", sep="")
    tempFile <- gsub(x=tempFile, pattern=" ", replacement="_")
    tempFile <- gsub(x=tempFile, pattern=":", replacement="-")
    if (is.null(ste)) {
        ste <- normalizePath(paste(tempFile, ".e", sep=""))
    }
    if (is.null(sto)) 
        sto <- normalizePath(paste(tempFile, ".o", sep=""))
    args1 <- paste(arguments, collapse=" ")
    catFile <- function(what, append=TRUE, fill=TRUE) {
        cat(what, file=tempFile, append=append, fill=fill)
    }
    catFile("#!/bin/bash", append=FALSE)
    if (!is.null(sourceBashrcBy)) 
        catFile(paste("source", sourceBashrcBy))
    catFile("")
    catFile(paste("#SBATCH -J", jobName))
    catFile(paste("#SBATCH -o", sto))
    catFile(paste("#SBATCH -e", ste))
    catFile("#SBATCH -n 1 ## total number of mpi tasks requested")
    catFile("#SBATCH -N 1 ## single node use")
    catFile(paste("#SBATCH -p", queue, "## queue"))
    catFile(paste("#SBATCH -t", time1))
    catFile("")
    catFile("echo \"Job started at:\"")
    catFile("echo $(date)")
    catFile(paste("Rscript ", sourceFile, " ", args1))
    catFile("echo \"Job ended at:\"")
    catFile("echo $(date)")
    command <- paste("sbatch", tempFile)
    if (doTalk) 
        print(command)
    if (wait < 1) 
        stop("Jobs may overwrite if wait<1 !")
    Sys.sleep(wait)
    while (TRUE) {
        jobsNum <- as.numeric(system("squeue -l -u $USER|wc -l", 
                                     intern=TRUE)) - 2
        if (jobsNum < maxJobs) 
            break
        Sys.sleep(5)
    }
    if (doSubmit) 
        system(command)
    if (doDeleteTempFile) 
        system(paste("rm -rf ", tempFile, sep=""))
    return(list(tempFile=tempFile, jobName=jobName, command=command))
}
