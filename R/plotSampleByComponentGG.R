#' plotSampleByComponentGG
#'
#' Plots a heatmap of the sample-by-component matrix. By default, feature
#' components are not clustered or ordered and listed in order provided and
#' scaling is applied across columns. Additional options can by provided to
#' alter the heat map rendered by [stats::heatmap()].
#'
#' @param object CNQuant or SigQuant class object containing a computed
#'   sample-by-component matrix
#' @param ... additional parameters passed to \link[stats]{heatmap}
#'
#' @return plot
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- createCNQuant(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- calculateFeatures(cnobj[1:10],method="drews")
#'   cnobj <- calculateSampleByComponentMatrix(cnobj)
#'   plotSampleByComponent(cnobj)
#' @seealso [getSampleByComponentGG()]
#' @export plotSampleByComponentGG
#'
plotSampleByComponentGG <- function(object=NULL,...){
    if(is.null(object)){
        stop("no object provided")
    }
    if(length(object@featFitting) == 0){
        stop("feature fitting not calculated")
    }

    plotData <- object@featFitting$sampleByComponent
    plotData <- apply(plotData,MARGIN = 2,FUN = function(x) (x - mean(x)) / sd(x))
    plotData <- as.data.frame(plotData) %>%
                    tibble::rownames_to_column("sample") %>%
                    tidyr::pivot_longer(cols = 2:ncol(.),values_to = "posterior",
                                        names_to = "component")

    ggplot2::ggplot(plotData) +
        ggplot2::geom_tile(ggplot2::aes(component,sample,fill=posterior)) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.border = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank())


    # stats::heatmap(x = plotData,
    #                scale="col",
    #                Colv = NA,
    #                ...)
}
