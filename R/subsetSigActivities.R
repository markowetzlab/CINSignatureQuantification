subsetSigActivities <- function(x,s){
    subSigActivities <- lapply(x, FUN = function(y){
        y <- y[which(rownames(y) %in% s),]
    })
    subSigActivities
}
