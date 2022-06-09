#' @rdname getSamplefeatures-methods
#' @aliases getSamplefeatures
setMethod("getSamplefeatures",signature = "CNQuant",function(object){
    object@samplefeatData
})
