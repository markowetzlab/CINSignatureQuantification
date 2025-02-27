#' plotDefinitions
#'
#' Plot the copy number feature components for a given sample in a SigQuant class object
#' containing computed sample by component. Default ordering by
#' decreasing exposure to signature CX1.
#'
#' @param object A SigQuant class object containing
#' @param cols Vector of colours of the same length as number of copy number
#'   components.
#' @param plot.dim vector length two giving the plotting grid dimensions
#' @return plot
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   t478 <- TCGA_478_Samples_SNP6_GOLD
#'   subsample <- unique(t478$sample)[1:10]
#'   t478 <- t478[t478$sample %in% subsample]
#'   cnobj <- quantifyCNSignatures(t478)
#'   plotDefinitions(cnobj)
#' @seealso [getDefinitions()]
#' @export plotDefinitions
#'
plotDefinitions <- function(object,cols=NULL,plot.dim=NULL,normalise=FALSE){
    if(is.null(object)){
        stop("No object provided, object should be a object of class SigQuant or SigQuant")
    }

    if(!inherits(object,"SigQuant")){
        stop("Object is not of class SigQuant or SigQuant")
    }

    if(nrow(object@backup.signatures) < 1){
        stop("No siganture by component matrix")
    }
    # check normalise arguement
    stopifnot(logical(normalise))
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
        n.feats <- sort(unique(gsub(x = colnames(object@featFitting$sampleByComponent),pattern = "[0-9]+$",replacement = "")))
        cols.per.feat <- unlist(lapply(n.feats,
                                       function(x){
                                           sum(grepl(pattern = x,
                                                     x = colnames(defs)))
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

    if(is.null(plot.dim)){
        n <- 1
        m <- 1
        side <- "n"
        while(prod(c(n,m)) < nrow(defs)){
            if(side == "n"){
                n <- n+1
                side <- "m"
            } else {
                m <- m+1
                side <- "n"
            }

        }

        setGrid <- c(n,m)
    } else {
        if(is.numeric(plot.dim)){
            if(length(plot.dim) == 2){
                if(prod(plot.dim) >= nrow(defs)){
                    setGrid <- plot.dim
                    n <- plot.dim[1]
                    m <- plot.dim[2]
                } else {
                    stop("plot.dim dimensions too small\n  plot.dim should be a c(n,m) vector where n x M >= number of signatures")
                }
            } else {
                stop("plot.dim incorrect length\n  plot.dim should be a c(n,m) vector where n x M >= number of signatures")
            }
        } else {
            stop("plot.dim not numeric\n  plot.dim should be a c(n,m) vector where n x M >= number of signatures")
        }
    }
    ## Include signature/feature normalisation used in paper figure
    if(normalise){
        ## First normalise per signature (to remove the comparable signature)
        theseNewSigs = apply(defs, 1, function(thisSig) {
            lSig = sapply(FEATS, function(thisFeat) {
                theseVals = thisSig[ grepl(thisFeat, names(thisSig)) ]
                theseNew = theseVals / sum(theseVals)
                # Catch edge case with only zeros and no weights (produce NaN)
                if(is.nan(sum(theseNew))) { 
                    altNew = rep(0, length(theseNew))
                    names(altNew) = names(theseNew)
                    theseNew = altNew
                }
                # Scale by numbers of components => Final sum of vector should be five (for five features)
                theseNew = theseNew * (length(theseNew)/NUMCOMP)
            } )
            vSig = unlist(lSig)
            return(vSig)
        } )
        defs = t(theseNewSigs)
        colnames(defs) = sapply(strsplit(colnames(defs), "\\."), function(x) x[[2]])
    }
    component <- feature <- value <- signature <- NULL
    tab <- as.data.frame(defs) %>%
        tibble::rownames_to_column(var = "signature") %>%
        tidyr::pivot_longer(cols = -1,names_to = "component") %>%
        dplyr::mutate(feature = gsub(pattern = "\\d+$",replacement = "",x = component)) %>%
        dplyr::mutate(signature = factor(x = signature,
                                         levels = stringr::str_sort(unique(signature),numeric = TRUE))) %>%
        dplyr::mutate(component = factor(x = component,
                                         levels = colnames(defs))) %>%
        dplyr::mutate(feature = factor(x = feature,
                                       levels = stringr::str_sort(unique(n.feats),numeric = TRUE)))


    ggplot2::ggplot(tab) +
        ggplot2::geom_col(ggplot2::aes(component,value,fill=feature)) +
        ggplot2::facet_wrap(. ~ signature,nrow = n,ncol = m) +
        ggplot2::scale_fill_manual(values = featcols) +
        ggplot2::scale_y_continuous(limits = c(0,1),expand = c(0,0)) +
        ggplot2::ylab(paste0("weight (",method,")")) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom",
                       axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.title.x = ggplot2::element_blank())
}
