# check chromosome input format and correct if needed
checkChromosomeFormat <- function(x=NULL){
    if(is.null(x)){
        stop("no chromosome col")
    }
    chrs <- c(seq.int(1:22),c("X"))
    # Catch chr prefix
    if(any(grepl(pattern = "chr",x = x))){
        message("checkChromosomeFormat: dropping 'chr' prefix")
        x <- gsub(pattern = "chr",replacement = "",x = x)
    }

    if(any(x %in% c("23","24"))){
        message("checkChromosomeFormat: formatting numerical sex chromsomes")
        x[x=="23"] <- "X"
        x[x=="24"] <- "Y"
    }

    if(any(!x %in% chrs)){
        message("checkChromosomeFormat: dropping unsupported chromsomes")
        x <- x[x %in% chrs]
    }
    # return
    x <- as.factor(x)
    return(x)
}
