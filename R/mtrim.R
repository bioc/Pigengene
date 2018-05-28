mtrim <- function (x){
    ## Removes white spaces from both sides of a string like a Refseq ID.
    gsub("^\\s+|\\s+$", "", x)
}
