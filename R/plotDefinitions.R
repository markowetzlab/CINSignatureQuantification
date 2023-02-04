#' plotDefinitions
#'
#' Plot the copy number feature components for a given sample in a SigQuant class object
#' containing computed sample by component. Default ordering by
#' decreasing exposure to signature CX1.
#'
#' @param object A SigQuant class object containing
#' @param cols Vector of colours of the same length as number of copy number
#'   components.
#' @return plot
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   t478 <- TCGA_478_Samples_SNP6_GOLD
#'   subsample <- unique(t478$sample)[1:10]
#'   t478 <- t478[t478$sample %in% subsample]
#'   cnobj <- quantifyCNSignatures(t478)
#'   activities <- plotDefinitions(cnobj)
#' @seealso [getSampleByComponent()]
#' @export plotDefinitions
#'
plotDefinitions <- function(object,cols=NULL){
    if(is.null(object)){
        stop("No object provided, object should be a object of class SigQuant or SigQuant")
    }

    if(!class(object) == c("SigQuant")){
        stop("Object is not of class SigQuant or SigQuant")
    }

    if(nrow(object@backup.signatures) < 1){
        stop("No sample by component matrix")
    }

    method <- object@signature.model
    defs <- object@backup.signatures

    if(!is.null(cols)){
        if(ncol(defs) != length(cols)){
            stop(paste0("Provided colour vector does not match number of components
                        \ncolours - ",length(cols),"; components - ",
                        ncol(defs)))
        }
    } else {
        n.components <- ncol(defs)
        n.feats <- names(object@featFitting$model)
        cols.per.feat <- unlist(lapply(n.feats,
                                       function(x){
                                           sum(grepl(pattern = x,
                                                     x = colnames(defs)))
                                       }))
        default.cols <- c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")
        switch(method,
               mac={
                   featcols <- default.cols
               },
               drews={
                   featcols <- default.cols[1:5]
               })
        cols <- rep(featcols,times=cols.per.feat)
        if(n.components != length(cols)){
            stop("Colour lengths do not match")
        }
    }

    par(mfrow = c(ceiling(sqrt(nrow(defs))),
                  floor(sqrt(nrow(defs)))))
    for(i in rownames(defs)){
    graphics::barplot(defs[i,],
                      main = paste0("CN signature ",i),
                      col = cols,
                      xlab = "component",
                      #names.arg=rep("",ncol(tabl)),
                      ylim = c(0,1),
                      ylab = paste0("weight (",method,")"),
                      axes=TRUE)
    }
    par(mfrow = c(1,1))
}
