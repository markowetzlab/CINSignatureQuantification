subsetfeatFitting <- function(x,s){
    sxc <- x$sampleByComponent
    sxc <- sxc[which(rownames(sxc) %in% s),]
    x$sampleByComponent <- sxc
    return(x)
}
