#' plotSampleByComponent
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
#' @seealso [getSampleByComponent()]
#' @export plotSampleByComponent
#'
plotSampleByComponent <- function(object=NULL,...){
    if(is.null(object)){
        stop("no object provided")
    }
    if(length(object@featFitting) == 0){
        stop("feature fitting not calculated")
    }
    plotData <- object@featFitting$sampleByComponent
    stats::heatmap(x = plotData,
                   scale="col",
                   Colv = NA,
                   ...)
}
