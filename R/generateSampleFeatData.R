generateSampleFeatData <- function(x){
    segCounts <- getSegCounts(x)
    ploidy <- getPloidyfeat(x)
    featData <- data.frame(
                        segCounts = segCounts,
                        ploidy = ploidy)
    return(featData)
}
