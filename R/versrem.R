versrem <- function (x) {
    ## Removes the version from (almost) Refseq IDs.
    y <- x
    if (length(unlist(strsplit(mtrim(x), "_"))) > 2) {
        y <- paste(unlist(strsplit(mtrim(x), "_"))[1:2], sep="_", 
                   collapse="_")
    }
    y
}
