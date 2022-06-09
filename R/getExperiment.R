#' @rdname getExperiment-methods
#' @aliases getExperiment
setMethod("getExperiment",signature = "CNQuant",function(object){
    object@ExpData
})
