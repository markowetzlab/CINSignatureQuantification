# Set input types for segtable columns
checkInputTypes <- function(x=NULL){
    if(is.null(x)){
        stop("no table")
    }
    x$chromsome <- as.character(x$chromosome)
    x$start <- as.numeric(x$start)
    x$end <- as.numeric(x$end)
    x$segVal <- as.numeric(x$segVal)
    x$sample <- as.factor(x$sample)
    return(x)
}
