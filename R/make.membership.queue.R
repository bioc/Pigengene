make.membership.queue <- function(moduleIDs, pigengene, modules){
    qu <- c()
    for(m1 in moduleIDs){
        qu <- c(qu, get.sorted.members(m1, pigengene, modules=modules ))}
    return(qu[order(rank(abs(qu)))])
}
