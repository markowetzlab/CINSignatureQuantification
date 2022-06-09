#' @rdname getSegments-methods
#' @aliases getSegments
setMethod("getSegments",signature = "CNQuant",function(object){
    segTable <- do.call(rbind,object@segments)
    data.frame(segTable,row.names = NULL)
})
