#' @rdname getSampleByComponent-methods
#' @aliases getSampleByComponent
setMethod("getSampleByComponent",signature = "CNQuant",function(object){
        object@featFitting$sampleByComponent
})
