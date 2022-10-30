#' @rdname getFeatures-methods
#' @aliases getFeatures
setMethod("getFeatures",signature = "CNQuant",function(object,feat=NULL){
    if(is.null(object@featData)){
        stop("no feature data available. Run `quantifyCNSignatures()` or
             `calculateFeatures()`")
    }
    if(is.null(feat)){
        object@featData
    } else {
        if(!all(feat %in% names(object@featData))){
            m <- feat[which(!feat %in% names(object@featData))]
            stop(paste0("Unknown feature name(s): ",paste0(m,collapse = ", "),"\n",
                 "available features: ",paste0(names(object@featData),collapse = ", ")))
        } else {
            object@featData[feat]
        }
    }
})
