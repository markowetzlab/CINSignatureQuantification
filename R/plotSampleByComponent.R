#' plotSampleByComponent
#'
#' Plots a heatmap of the sample-by-component matrix
#'
#' @param object a CNQuant or SigQuant class object
#' @param ... additional parameters passed to \link[stats]{heatmap}
#'
#' @return plot
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
    stats::heatmap(x = plotData,...)
}
