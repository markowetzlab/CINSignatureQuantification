subsetCNfeatures <- function(x,s){
    subCNfeats <- lapply(x, FUN = function(y){
        y <- y[y$ID %in% s,]
    })
    subCNfeats
}
