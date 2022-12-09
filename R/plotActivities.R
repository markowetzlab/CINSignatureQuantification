#' plotActivities
#'
#' Plot the copy number signature activities for a given SigQuant class object
#' containing copy number signature activities/exposures. Default ordering by
#' decreasing exposure to signature CX1.
#'
#' @param object A SigQuant class object containing
#' @param type type of copy number signature matrix to return 'raw', 'norm',
#'   'threshold', or 'scaled' (default: threshold).
#' @param cols Vector of colours of the same length as number of copy number
#'   signatures.
#' @return plot
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   t478 <- TCGA_478_Samples_SNP6_GOLD
#'   subsample <- unique(t478$sample)[1:10]
#'   t478 <- t478[t478$sample %in% subsample]
#'   cnobj <- quantifyCNSignatures(t478)
#'   activities <- plotActivities(cnobj,type="threshold")
#' @seealso [getActivities()]
#' @export plotActivities
#'
plotActivities <- function(object=NULL,type="threshold",cols=NULL){
    if(is.null(object)){
        stop("No object provided")
    }
    if(!class(object) == "SigQuant"){
        stop("Object is not of class SigQuant")
    }

    method <- object@signature.model
    types <- c("raw","norm","threshold","scaled")
    if(!type %in% c("raw","norm","threshold","scaled")){
        stop(paste0("Unknown signature type provided. Provide either ",paste0(types,collapse = ", ")))
    }
    switch(method,
           mac={
               if(!type %in% c("raw","threshold")){
                   stop("Signature method only supports raw or normalised threshold-adjusted")
               }
               switch(type,
                      raw={
                          plotdata <- object@activities$rawAct0
                      },
                      threshold={
                          plotdata <- object@activities$thresholdAct2
                      })
           },
           drews={
               switch(type,
                      raw={
                          plotdata <- object@activities$rawAct0
                      },
                      norm={
                          plotdata <- object@activities$normAct1
                      },
                      threshold={
                          plotdata <- object@activities$thresholdAct2
                      },
                      scaled={
                          plotdata <- object@activities$scaledAct3
                      })
           })
    if(!is.null(cols)){
        if(ncol(plotdata) != length(cols)){
            stop(paste0("Provided colour vector does not match number of signatures.
                        \ncolours - ",length(cols),"; signatures - ",ncol(plotdata)))
        }
    } else {
        switch(method,
               mac={
                   cols <- c(s1="#1B9E77",s2="#D95F02",s3="#7570B3",
                             s4="#E7298A",s5="#66A61E",s6="#E6AB02",
                             s7="#A6761D")
                   cols <- cols[colnames(plotdata)]
               },
               drews={
                   cols <- c(CX1='#a6cee3',CX2='#1f78b4',CX3='#b2df8a',
                             CX4='#33a02c',CX5='#fb9a99',CX6='#e31a1c',
                             CX7='#fdbf6f',CX8='#ff7f00',CX9='#cab2d6',
                             CX10='#6a3d9a',CX11='#ffff99',CX12='#b15928',
                             CX13='#8dd3c7',CX14='#ffffb3',CX15='#bebada',
                             CX16='#fb8072',CX17='#80b1d3')
                   cols <- cols[colnames(plotdata)]
               })

    }

    clms <- 1
    l.pos <- -0.1

    if(nrow(plotdata) > 1){
        plotdata <- plotdata[order(plotdata[,1],decreasing = T),]
    }
    tabl <- t(as.matrix(plotdata))

    graphics::par(mar=c(5, 4, 4, 8), xpd=TRUE)
    graphics::barplot(tabl,
            main = paste0("Signature activities (","method: ",method,")"),
            col = cols,
            xlab = "sample",
            names.arg=rep("",ncol(tabl)),
            ylab = paste0("activity (",type,")"),
            axes=TRUE)
    graphics::legend(x= "topright",
           inset=c(l.pos, 0),
           legend = rownames(tabl),
           fill=cols,
           cex=0.7,
           ncol=clms,
           x.intersp = 0.5,
           y.intersp = 0.8,
           box.col=NA)
}
