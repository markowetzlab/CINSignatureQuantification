getSegCounts <- function(x){
    segCounts <- unlist(lapply(x,nrow))
    return(segCounts)
}
