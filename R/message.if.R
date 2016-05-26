message.if <- function(me=NULL, verbose=0){
    ## Prints the message(s) if verbose is more than 0.
    if(verbose>0)
        for(m1 in me)
            message(m1)
}
