gene.mapping <- function (
    ids, inputType="REFSEQ", outputType="SYMBOL", leaveNA=TRUE, 
    inputDb="Human", outputDb=inputDb, verbose=0)
{
    ## inputDb: Input database.
    ##^It can be org.Hs.eg.db for human and org.Mm.eg.db for mouse.
    ##
    addAt <- FALSE
    if (outputType=="ENTREZIDat") {
        addAt <- TRUE
        outputType <- "ENTREZID"
    }
    if (inputType=="Auto") {
        if (length(grep(ids[1], pattern="NM_|NR_"))) 
            inputType <- "REFSEQ"
        if (inputType=="Auto") 
            stop("inputType could not be determined automatically!")
    }
    if (inputType=="ENTREZIDat") {
        ids <- gsub(ids, pattern="_at", replacement="")
        inputType <- "ENTREZID"
    }
    key <- mtrim(ids)
    ##Remove the version number, i.e. anything after the second "_". 
    if(inputType=="REFSEQ")
        key <- unlist(lapply(X=key, FUN=versrem))
    key0 <- key
    ## The default db: ## AFTER ACCEPTANCE
    if(class(inputDb)=="character"){
        if(inputDb=="Human")
            inDbChar <- "org.Hs.eg.db"
    } else {
        inDbChar <- inputDb$packageName
    }
    if(class(outputDb)=="character"){
        if(outputDb=="Human")
            outDbChar <- "org.Hs.eg.db"
    } else {
        outDbChar <- outputDb$packageName
    }
    inputTypeOrig <- inputType
    ## Check the availability of the database:
    if("org.Hs.eg.db" %in% c(inDbChar,outDbChar) & !require(org.Hs.eg.db))
       stop("org.Hs.eg.db package is required by gene.mapping()!")
    if("org.Mm.eg.db" %in% c(inDbChar,outDbChar) & !require(org.Mm.eg.db))
       stop("org.Mm.eg.db package is required by gene.mapping()!")
    ## Change char to database package:
    inputDb <- get(inDbChar)
    outputDb <- get(outDbChar)
    ## Mapping between species:
    if(inDbChar!= outDbChar){
        ens <- gene.mapping(ids=key, inputType=inputType, outputType="ENSEMBL", 
                            inputDb=inputDb, verbose=verbose-1, leaveNA=TRUE)[, 3]
        ends <- ens[!is.na(ens)]
        inputType <- "ENSEMBL"
        hm <- hu.mouse(verbose=verbose)
        if(inDbChar=="org.Hs.eg.db" & outDbChar=="org.Mm.eg.db")
            key <- hm$hu2mouse[ens]
        if(inDbChar=="org.Mm.eg.db" & outDbChar=="org.Hs.eg.db")
            key <- hm$mouse2hu[ens]
    } ## else: input and output databases are the same. 
    ##
    if(require(AnnotationDbi)){
        q1B <- AnnotationDbi::select(outputDb, keys=key, columns=outputType,
                                     keytype=inputType)
    } else {
        stop("biomaRt package is required by gene.mapping()!")
    }
    colnames(q1B)[1] <- inputType ## important when inputType=outputType.
    output2 <- as.character(q1B[match(key, q1B[, inputType]), outputType])
    f1 <- cbind(key0, output1=output2, output2=output2)
    colnames(f1)[1] <- "input"
    if (!leaveNA) {
        inds <- which(is.na(f1[, "output2"]))
        nms <- as.character(ids[inds])
        if (length(inds) > 0) {
            f1[inds, "output2"] <- nms
        }
    }
    if (addAt) 
        f1[, "output2"] <- paste(f1[, "output2"], "_at", sep="")
    rownames(f1) <- ids
    return(f1)
}
