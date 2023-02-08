#' plotSpectrum
#'
#' Plot the copy number feature components for a given sample in a CNQuant class object
#' containing computed sample by component. Default ordering by
#' decreasing exposure to signature CX1.
#'
#' @param object A CNQuant class object containing
#' @param sample vector of length 1 containing either a sample name or sample
#'   index
#' @param cols Vector of colours of the same length as number of copy number
#'   components.
#' @return plot
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   t478 <- TCGA_478_Samples_SNP6_GOLD
#'   subsample <- unique(t478$sample)[1:10]
#'   t478 <- t478[t478$sample %in% subsample]
#'   cnobj <- quantifyCNSignatures(t478)
#'   plotSpectrum(cnobj,sample = 1)
#' @seealso [getSampleByComponent()]
#' @export plotSpectrum
#'
plotSpectrum <- function(object,sample=NULL,cols=NULL){
    if(is.null(object)){
        stop("No object provided, object should be a object of class CNQuant or SigQuant")
    }

    if(!class(object) %in% c("CNQuant","SigQuant")){
        stop("Object is not of class CNQuant or SigQuant")
    }

    if(length(object@featFitting) < 1){
        stop("No sample by component matrix")
    }

    if(is.null(sample)){
        stop("No sample specified; sample should be an integer index or name of sample contained within the provided object")
    }

    if(!is.numeric(sample) & !is.character(sample)){
        stop("Unknow sample value provided; sample should be an integer index or name of sample contained within the provided object")
    }

    samp <- getSamples(object = object)
    samp.len <- length(samp)
    if(is.numeric(sample)){
        if(sample > samp.len){
            stop(paste0("Sample index is out of bounds; Object contains ",samp.len," samples"))
        }
    }

    if(is.character(sample)){
        if(!sample %in% samp){
            stop("Sample was not found in the object provided")
        }
    }

    method <- object@featFitting$method
    sxc <- object@featFitting$sampleByComponent

    if(!is.null(cols)){
        if(ncol(sxc) != length(cols)){
            stop(paste0("Provided colour vector does not match number of components
                        \ncolours - ",length(cols),"; components - ",
                        ncol(sxc)))
        }
    } else {
        n.components <- ncol(sxc)
        n.feats <- names(object@featFitting$model)
        cols.per.feat <- unlist(lapply(n.feats,
                                       function(x){
                                           sum(grepl(pattern = x,
                                                     x = colnames(sxc)))
                                       }))
        default.cols <- c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")
        names(default.cols) <- c("segsize","changepoint","bp10MB","bpchrarm","osCN","copynumber")
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

    samp.name <- ifelse(is.numeric(sample),samp[sample],sample)
    object <- object[samp.name]
    tabl <- object@featFitting$sampleByComponent
    tabl <- tabl / sum(tabl)

    # graphics::barplot(tabl,
    #                   main = paste0("sample-by-component spectrum (",samp.name,")"),
    #                   col = cols,
    #                   xlab = "component",
    #                   #names.arg=rep("",ncol(tabl)),
    #                   ylab = paste0("Sum-of-posterior (",method,")"),
    #                   axes=TRUE)
    component <- feature <- value <- NULL
    tab <- as.data.frame(tabl) %>%
        tibble::rownames_to_column(var = "component") %>%
        dplyr::mutate(feature = gsub(pattern = "\\d+$",replacement = "",x = component)) %>%
        dplyr::rename("value"="tabl") %>%
        dplyr::mutate(component = factor(x = component,
                                         levels = colnames(sxc))) %>%
        dplyr::mutate(feature = factor(x = feature,
                                       levels = stringr::str_sort(unique(n.feats),numeric = TRUE)))


    ggplot2::ggplot(tab) +
        ggplot2::geom_col(ggplot2::aes(component,value,fill=feature)) +
        ggplot2::scale_fill_manual(values = featcols) +
        ggplot2::scale_y_continuous(limits = c(0,1),expand = c(0,0)) +
        ggplot2::ylab(paste0("weight (",method,")")) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom",
                       axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.title.x = ggplot2::element_blank())
}
