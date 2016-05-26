get.module.data.file <- function(
    resultPath, moduleNum, use.Hartemink, partition=NULL) 
{
    result <- list()
    modulePath <- file.path(resultPath, paste("mod", moduleNum, sep=""))
    if (!is.null(partition))
        modulePath <- file.path(modulePath, paste("part", partition, sep=""))
    dir.create(path=modulePath, recursive=TRUE, showWarnings=FALSE)
    moduleFile <-  file.path(modulePath, paste("module", moduleNum, "H", 
                                               as.numeric(use.Hartemink), ".RData", sep=""))
    result[["modulePath"]] <- modulePath
    result[["moduleFile"]] <- moduleFile
    return(result)
}
