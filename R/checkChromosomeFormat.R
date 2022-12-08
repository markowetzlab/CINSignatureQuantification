# check chromosome input format and correct if needed
checkChromosomeFormat <- function(x=NULL){
    if(is.null(x)){
        stop("no chromosome col")
    }

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
    return(x)
}
