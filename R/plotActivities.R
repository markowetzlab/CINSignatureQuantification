#' plotActivities
#'
#' Plot the copy number signature activites for a given CNQuant or SigQuant class object
#' containing copy number signature activities/exposures. Default ordering by
#' signature CX1
#'
#' @param object A SigQuant class object
#'
#' @return plot
#' @export plotActivities
#'
plotActivities <- function(object){
    if(is.null(object)){
        stop("No object provided, object should be a object of class CNQuant or SigQuant")
    }
    if(!class(object) == "SigQuant"){
        stop("Object is not of class SigQuant")
    }
    ## May need to change which matrix is used
    if(object@signature.model == "drews"){
        plotdata <- object@activities$thresholdAct2
        cols <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
                 '#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
                 '#cab2d6','#6a3d9a','#ffff99','#b15928',
                 '#8dd3c7','#ffffb3','#bebada','#fb8072',
                 '#80b1d3')
        clms <- 1
        l.pos <- -0.1
    } else {
        plotdata <- object@activities$thresholdAct2
        cols <- c("#1B9E77","#D95F02","#7570B3","#E7298A",
                  "#66A61E","#E6AB02","#A6761D")
        clms <- 1
        l.pos <- -0.1
    }

    plotdata <- plotdata[order(plotdata[,1],decreasing = T),]
    tabl <- t(as.matrix(plotdata))

    par(mar=c(5, 4, 4, 8), xpd=TRUE)
    barplot(tabl,
            main = paste0("Signature activities (","method: ",object@signature.model,")"),
            col = cols,
            xlab = "sample",
            names.arg=rep("",ncol(tabl)),
            ylab = "relative exposure",
            axes=TRUE)
    legend(x= "topright",
           inset=c(l.pos, 0),
           legend = rownames(tabl),
           fill=cols,
           cex=0.7,
           ncol=clms,
           x.intersp = 0.5,
           y.intersp = 0.8,
           box.col=NA)
}
