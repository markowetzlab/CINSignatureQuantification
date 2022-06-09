getPloidyfeat <- function(x){
    featploidy <- unlist(lapply(x,FUN = function(y){
        segLen<-(as.numeric(y$end)-as.numeric(y$start))
        ploidy<-sum((segLen/sum(segLen))*as.numeric(y$segVal))
    }))
    featploidy <- round(featploidy,digits = 3)
    return(featploidy)
}
