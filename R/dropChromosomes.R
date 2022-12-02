dropChromosomes <- function(x=NULL){
    if(is.null(x)){
        stop("no seg table")
    }
    chrs <- c(seq.int(1:22),c("X"))
    if(any(!x$chromosome %in% chrs)){
        message("checkChromosomeFormat: dropping unsupported chromsomes")
        x <- x[x$chromosome %in% chrs,]
    }
    x$chromosome <- factor(x$chromosome,levels=chrs)
    return(x)
}
