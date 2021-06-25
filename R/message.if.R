message.if <- function(me=NULL, verbose=0, txtFile=NULL, append=TRUE, ...){
    ## Prints the message(s) if verbose is more than 0.
    if(verbose>0){
        for(m1 in me){
            message(m1)
        }
    }
    if(!is.null(txtFile)){
        capture.output(me, file=txtFile, append=append, ...) 
    }
}
