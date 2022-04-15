hu.mouse <- function(host="www.ensembl.org", verbose=0, mouseHomologFilter="with_mmusculus_homolog"){
    ##mouseHomologFilter used to be ""with_homolog_mmus"" in older versions of biomaRt.
    ##
    ## This function uses Biomart Ensembel to compute a table with human ensembl genes
    ##^ at the first column and their homologs at the second column.
    result <- list()
    ##
    message.if("Using Biomart Ensembel to map homolog genes...", verbose=verbose)
    attributes <- c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene")
    if (require(biomaRt)){
        ensembl <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                                    verbose=max(verbose-2, 0), 
                                    dataset="hsapiens_gene_ensembl", host=host)
        homologs <- biomaRt::getBM(attributes=attributes, values=TRUE,
                                   filters=mouseHomologFilter, uniqueRows=TRUE, 
                                   mart=ensembl, verbose=max(verbose-2, 0))
    } else {
        stop("biomaRt package is required by hu.mouse()!")
    }
    message.if("Biomart Ensembel DONE.", verbose=verbose)
    expected <- "ensembl_gene_id-mmusculus_homolog_ensembl_gene"
    if(paste(colnames(homologs), collapse="-") != expected)
        stop("Unexpected output from biomaRt::getBM !")
    colnames(homologs) <- c("human", "mouse")
    ## Mouse to human:
    ids <- unique(homologs[, "mouse"])
    mouse2huMatrix <- homologs[match(ids, homologs[, "mouse"]), ]
    mouse2hu <- mouse2huMatrix[, "human"]
    names(mouse2hu) <- mouse2huMatrix[, "mouse"]
    ## Human to mouse:
    ids <- unique(homologs[, "human"])
    hu2mouseMatrix <- homologs[match(ids, homologs[, "human"]), ]
    hu2mouse <- hu2mouseMatrix[, "mouse"]
    names(hu2mouse) <- hu2mouseMatrix[, "human"]
    multiHu <- nrow(homologs)-length(hu2mouse)
    multiMouse <- nrow(homologs)-length(mouse2hu)
    message.if(paste("The number of non-unique human-mouse genes=", multiHu), 
               verbose=verbose)
    message.if(paste("The number of non-unique mouse-human genes=", multiMouse), 
               verbose=verbose)
    warning("Some genes with multiple homologs are selected arbitrarily.")
    ##
    result[["homologs"]] <- homologs
    result[["hu2mouse"]] <- hu2mouse
    result[["mouse2hu"]] <- mouse2hu
    result[["multiMouse"]] <- multiMouse
    result[["multiHu"]] <- multiHu
    return(result)
}
