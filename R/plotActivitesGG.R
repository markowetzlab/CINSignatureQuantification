#' plotActivitiesGG
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
#' @export plotActivitiesGG
#'
plotActivitiesGG <- function(object,type="threshold",cols=NULL){
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

    ## allow for single sample plotting
    if(ncol(as.data.frame(plotdata)) == 1){
        plotdata <- t(plotdata)
        rownames(plotdata) <- "sample"
    }

    plotdata <- as.data.frame(plotdata) %>%
                    tibble::rownames_to_column(var = "sample") %>%
                    tidyr::pivot_longer(cols = 2:ncol(.),
                                        names_to = "signature",
                                        values_to = "activity") %>%
                    dplyr::mutate(signature = factor(x = signature,
                                                     levels = paste0(names(cols)))) %>%
                    dplyr::arrange(signature,desc(activity)) %>%
                    dplyr::mutate(sample = factor(x = sample,levels = unique(sample)))

    sigPlot <- ggplot2::ggplot(plotdata) +
                ggplot2::geom_col(ggplot2::aes(x = sample,y = activity,fill=signature),
                                  position = ggplot2::position_fill(reverse = TRUE),color="grey20") +
                ggplot2::scale_fill_manual(values = cols) +
                ggplot2::scale_y_continuous(expand = c(0,0)) +
                ggplot2::scale_x_discrete(expand = c(0,1)) +
                ggplot2::ggtitle(label = paste0("Signature activities (","method: ",method,")")) +
                ggplot2::ylab(label = paste0("activity (",type,")")) +
                ggplot2::theme_bw() +
                ggplot2::theme(legend.position = "right",
                               panel.border = ggplot2::element_blank(),
                               panel.grid = ggplot2::element_blank(),
                               axis.text.x = ggplot2::element_blank(),
                               axis.ticks.x = ggplot2::element_blank(),
                               axis.line.y.left = ggplot2::element_line())
    sigPlot
}
