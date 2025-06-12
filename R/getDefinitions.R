#' @rdname getDefinitions-methods
#' @aliases getDefinitions
setMethod("getDefinitions",signature = "SigQuant",function(object,normalise=FALSE){
    # check normalise arguement
    stopifnot(is.logical(normalise))

    method <- object@signature.model
    defs <- object@backup.signatures

    if(normalise){
        feats <- names(object@featData)
        if(method %in% c("drews")){
            feats <- feats[1:5]
        }
        ## First normalise per signature (to remove the comparable signature)
        theseNewSigs = apply(defs, 1, function(x) {
            lSig = sapply(feats, function(y) {
                theseVals = x[ grepl(y, names(x)) ]
                theseNew = theseVals / sum(theseVals)
                # Catch edge case with only zeros and no weights (produce NaN)
                if(is.nan(sum(theseNew))) {
                    altNew = rep(0, length(theseNew))
                    names(altNew) = names(theseNew)
                    theseNew = altNew
                }
                # Scale by numbers of components => Final sum of vector should be five (for five features)
                theseNew = theseNew * (length(theseNew)/ncol(defs))
            } )
            vSig = unlist(lSig)
            return(vSig)
        } )
        defs = t(theseNewSigs)
        colnames(defs) = sapply(strsplit(colnames(defs), "\\."), function(x) x[[2]])
    }

    return(defs)
})
