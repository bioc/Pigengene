get.enriched.pw <- function(genes, idType, pathwayDb, ont=c("BP", "MF", "CC"),
                            Org="Human", OrgDb=NULL, outPath, 
                            pvalueCutoff=0.05, pAdjustMethod="BH", fontSize=14, 
                            verbose=0){
    ## Isha wrote this function to perform functional enrichment analysis of selected
    ## genes.
    ## Input:
    ## genes: A character vector of genes for which pathway over representation analysis
    ## to be done.
    ## idType: A character string describing the type of input gene ID e.g.,
    ## "ENTREZID", "REFSEQ", "SYMBOL"
    ## pathwayDb: A character vector of enrichment database to be used e.g.,
    ## "GO", "KEGG", "REACTOME", or "NCG".
    ## ont: GO ontology terms to be analysed e.g. "BP", "MF" or "CC". Default is
    ## all three.
    ## Org: A character string equal to "Human" or "Mouse" determining the reference
    ## organism to be used.
    ## OrgDb: The reference data base to be used. Use \code{org.Ce.eg.db} for
    ## 'Celegans' org. For default org 'Human' and for org 'Mouse' the
    ## databses will be set to \code{org.Hg.eg.db} and \code{org.Mm.eg.db}
        ## respectively.}
    ## outPath: A file path where results will be saved.

    ## This function should return at a minimum, noEnrichment and the list of enrichment results. H?

    message.if(me="Getting enriched pathways...", verbose=verbose)
    ## QC:
    if(!is.null(Org)){
        if(!is.null(OrgDb)){
            stop("Provide a value for either 'OrgDb' or 'Org', but not both. 
             The other value will be updated based on provided value! Exactly
             one of the two must be NULL.")
        }
        ## Org is not NULL, OrgDb=NULL
        ## Setting values for Org and OrgDb
        if(any(Org %in% c("Human", "Mouse")) && is.null(OrgDb)){
            OrgDb <- get(type2char(type1=Org))
        } else{
            stop("Org cannot take any values other than 'Human' or 'Mouse'!")
        }       
    } ## OrgDb is not NULL any more
    
    loaded <- try(require(OrgDb$packageName, character.only=TRUE))
    if(inherits(loaded, "try-error") | !loaded)
       stop("The package determined by OrgDb input is not installed!")
    if(!any(ont %in% c("BP", "MF", "CC")))
        stop("ont cannot take values other than BP, MF or CC!")
    if(inherits(OrgDb, "character"))
        stop("OrgDb cannot be a vector!")

    ## Save results
    result <- list()
    resultPath <- file.path(outPath, "enrichResults")
    
    if(!any(Org %in% c("Human", "Mouse")) && !is.null(OrgDb)){
        Org <- as.character(DBI::dbGetQuery(AnnotationDbi::dbconn(OrgDb),
                                            "SELECT value FROM metadata WHERE name='SPECIES'"))
    }

    ## Identifying KEGG organism input
    keggOrg <- "hsa"
    if(Org=="Mouse"){
        keggOrg <- "mmu"
        ## refOut <- list("Human", "Mouse")
    }

    if(inherits(genes, "list")){
        ## QC
        if(length(names(genes))==0)
            stop("If genes is a list, it must have names!")
        resultS <- list()
        for(name1 in names(genes)){
            message.if(me=paste("Analysing: ", name1), verbose=verbose-1)
            genes1 <- genes[[name1]]
            resultPath <- file.path(outPath, paste0(name1, "_enrichResults"))
            resultS[[name1]] <- get.enriched.pw(genes=genes1, idType=idType,
                                                 pathwayDb=pathwayDb,
                                                 ont=ont, Org=NULL, OrgDb=OrgDb,
                                                 outPath=resultPath,
                                                 pvalueCutoff=pvalueCutoff,
                                                 pAdjustMethod=pAdjustMethod,
                                                 verbose=verbose)
        }
        return(resultS)
    }

    dir.create(resultPath, showWarnings=FALSE,  recursive=TRUE)
 
    if(idType!="ENTREZID"){
        message.if(me="Mapping input gene IDs to ENTREZID", verbose=verbose-1)

        ## Converting genes inputID type
        targetConv <- gene.mapping(ids=na.omit(as.character(genes)),
                                   inputType=idType,
                                   outputType="ENTREZID",
                                   inputDb=OrgDb) ##outputDb=Org
        targetConv <- as.data.frame(targetConv, stringsAsFactors=FALSE)
        colnames(targetConv) <- c(idType, "ENTREZID")
        genes <- na.omit(targetConv$ENTREZID)
    }
    genes <- na.omit(genes)

    tables <- list()
    enrichmentS <- list()
    for(e1 in pathwayDb){
        if(e1=="GO"){
            for(g1 in ont){
                enriched <- enrichGO(gene=genes, OrgDb=OrgDb, keyType="ENTREZID",
                                     ont=g1, pvalueCutoff=pvalueCutoff,
                                     pAdjustMethod=pAdjustMethod)
                enrichedSimple <- simplify(enriched)
                ##^ all child branches of a pathway are represented by one parent pathway
                enrichmentS[[paste0(e1,"-",g1)]] <- enrichedSimple
                tables[[paste0(e1,"-",g1)]] <- data.frame(enrichedSimple)
            }
        }

        if(e1=="KEGG"){
            eKEGG <- enrichKEGG(gene=genes, organism=keggOrg,
                                pvalueCutoff=pvalueCutoff,
                                pAdjustMethod=pAdjustMethod)
            enrichmentS[[e1]] <- eKEGG
            tables[[e1]] <- data.frame(eKEGG)
        }

        if(e1=="REACTOME"){
            eREACT <- enrichPathway(gene=genes, organism=tolower(Org),
                                    pvalueCutoff=pvalueCutoff,
                                    pAdjustMethod=pAdjustMethod, readable=TRUE)
            enrichmentS[[e1]] <- eREACT
            tables[[e1]] <- data.frame(eREACT)
        }

        if(e1=="NCG"){
            eNCG <- enrichNCG(gene=genes, pvalueCutoff=pvalueCutoff,
                              pAdjustMethod=pAdjustMethod, readable=TRUE)
            enrichmentS[[e1]] <- eNCG
            tables[[e1]] <- data.frame(eNCG)
        }
    }
    
    ## Adding gene symbols
    for (tableInd in 1:length(tables)) {
        table1 <- tables[[tableInd]]
        if (nrow(table1)==0)
            next
        geneSymbols <- c()
        for (ri in 1:nrow(table1)) {
            overlapGenes <- table1[ri , "geneID"]
            geneIds <- unlist(strsplit(overlapGenes, split="/"))
            mapped <- gene.mapping(ids=geneIds, inputDb="Human" ,
                                   inputType="ENTREZID", outputType="SYMBOL", verbose=verbose-2)
            mapped <- mapped[ ,"output2"]
            mapped <- mapped[!is.na(mapped)]
            geneSymbols <- c(geneSymbols, paste(mapped, collapse=" ")) 
        }
        table1 <- cbind(table1, geneSymbols)
        ## Moving the "Count" column to the last column
        table1 <- cbind(table1[ ,colnames(table1)!= "Count"], Count=table1[ ,"Count"])
        tables[[tableInd]] <- table1
    }

    ## Excluding empty databases
    subListData <- Filter(function(x) nrow(x) > 0, tables)
    subEnrichedList <- enrichmentS[names(subListData)]
    noEnrichment <- names(tables)[!names(tables) %in%
                                  names(subListData)]
    if(length(noEnrichment) != 0)
        message.if(me=paste("No enriched pathways found using", noEnrichment, "database"), 
                   verbose=verbose)
    result[["EnrichResults"]] <- subEnrichedList
    result[["noEnrichment"]] <- noEnrichment

    ## Save the pathway data in excel sheet
    xlsFile <- file.path(resultPath, paste0("ORA_results_", pvalueCutoff, "_",
                                            pAdjustMethod, ".xlsx"))
    if(length(subListData) > 0){
        if(file.exists(xlsFile)){
            file.remove(xlsFile)
        }
        write.xlsx(x=subListData, file=xlsFile)
    }

    ## Plot
    for(l1 in names(subEnrichedList)){
        dotplot(subEnrichedList[[l1]], showCategory=30, orderBy="GeneRatio") +
                ggtitle(paste0(l1, "-Enriched Pathways")) +
                theme(plot.title=element_text(size=fontSize, face="bold", hjust=0.5),
                axis.text.y=element_text(size=fontSize, face="bold"),
                axis.text.x=element_text(size=12, face="bold"))
        plotFile <- file.path(resultPath, paste0(l1, "_ORA_", pvalueCutoff, "_",
                                                 pAdjustMethod, ".png"))
        ggsave(plotFile, device="png", width=10, height=10, dpi="retina")
    }

    message.if(me=paste("Plots and an excel file are saved at:", resultPath), 
               verbose=verbose)


    return(result)
}
