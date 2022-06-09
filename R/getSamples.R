#' @rdname getSamples-methods
#' @aliases getSamples
setMethod("getSamples",signature = "CNQuant",function(object){
    rownames(object@samplefeatData)
})
