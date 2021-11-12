gene.mapping <- function(ids, inputType="REFSEQ", outputType="SYMBOL", leaveNA=TRUE, 
                         inputDb="Human", outputDb=inputDb, verbose=0){
    ## inputDb: Input database.
    ##^It can be 'Human' or org.Hs.eg.db for human and 'Mouse' or org.Mm.eg.db for mouse.
 
    ## Sinlge or multiple output types?
    if(length(outputType) >1 | length(outputDb) >1){
        res <- c()
        ## Creating a column for each desired output type or DB:
        if(inherits(outputDb, "list")){
            outputDbList <- outputDb
        } else { ## It is a single Db,
            outputDbList <- list(outputDb)
        }
        for(od in outputDbList){
            for(ot in outputType){
                odChar <- type2char(type1=od)
                mapTo <- paste(odChar, ot, sep="-")
                message.if(cat("Mapping to: ", mapTo, "\n"), verbose=verbose)
                mapped <- gene.mapping(ids=ids, inputType=inputType, outputType=ot, leaveNA=TRUE, 
                                       inputDb=inputDb, outputDb=od, verbose=verbose-1)
                mapped <- mapped[,"output1", drop=FALSE]
                res <- cbind(res, mapped)
                colnames(res)[ncol(res)] <- mapTo
            }
        }
        return(res)
    }

    ## Cleaning:
    addAt <- FALSE
    if(outputType=="ENTREZIDat") {
        addAt <- TRUE
        outputType <- "ENTREZID"
    }
    if(inputType=="Auto") {
        if(length(grep(ids[1], pattern="NM_|NR_"))) 
            inputType <- "REFSEQ"
        if(inputType=="Auto") 
            stop("inputType could not be determined automatically!")
    }
    if(inputType=="ENTREZIDat") {
        ids <- gsub(ids, pattern="_at", replacement="")
        inputType <- "ENTREZID"
    }
    key <- mtrim(ids)
    ##Remove the version number, i.e. anything after the second "_". 
    if(inputType=="REFSEQ")
        key <- unlist(lapply(X=key, FUN=versrem))
    key0 <- key
    ## The default db: ## AFTER ACCEPTANCE
    inDbChar <- type2char(type1=inputDb)
    outDbChar <- type2char(type1=outputDb)

    inputTypeOrig <- inputType
    ## Check the availability of the database:
    if("org.Hs.eg.db" %in% c(inDbChar, outDbChar) & !require(org.Hs.eg.db))
       stop("org.Hs.eg.db package is required by gene.mapping()!")
    if("org.Mm.eg.db" %in% c(inDbChar, outDbChar) & !require(org.Mm.eg.db))
       stop("org.Mm.eg.db package is required by gene.mapping()!")
    ## Change char to database package:
    inputDb <- get(inDbChar)
    outputDb <- get(outDbChar)

    ## QC:
    possibleKeys <- AnnotationDbi::keytypes(inputDb)
    possibleCols <- AnnotationDbi::columns(outputDb)
    if(!inputType %in% possibleKeys)
        stop(paste(inputType, "not in possible keys:", paste(possibleKeys, collapse=", ")))
    if(!outputType %in% possibleCols)
        stop(paste(outputType, "not in possible cols:", paste(possibleKeys, collapse=", ")))
    
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
    if(!leaveNA) {
        inds <- which(is.na(f1[, "output2"]))
        nms <- as.character(ids[inds])
        if(length(inds) > 0) {
            f1[inds, "output2"] <- nms
        }
    }
    if(addAt) 
        f1[, "output2"] <- paste(f1[, "output2"], "_at", sep="")
    rownames(f1) <- ids
    return(f1)
}
