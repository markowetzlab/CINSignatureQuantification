#' @rdname getSamplefeatures-methods
#' @aliases getSamplefeatures
setMethod("getSamplefeatures",signature = "CNQuant",function(object){
    as.data.frame(object@samplefeatData)
})
