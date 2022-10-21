#' @rdname getSignatures-methods
#' @aliases getSignatures
setMethod("getSignatures",signature = "SigQuant",function(object,type="threshold"){
    types <- c("raw","norm","threshold","scaled")
    method <- getExperiment(object)@feature.method
    if(!type %in% c("raw","norm","threshold","scaled")){
        stop(paste0("Unknown signature type provided. Provide either ",paste0(types,collapse = ", ")))
    }
    switch(method,
            mac={
               if(!type %in% c("raw","threshold")){
                   stop("Signature method only supports raw or normalised threshold-adjusted")
                }
                switch(type,
                       raw={
                           return(as.matrix(object@activities$rawAct0))
                       },
                       threshold={
                           return(as.matrix(object@activities$thresholdAct2))
                       })
           },
           drews={
               switch(type,
                      raw={
                          return(as.matrix(object@activities$rawAct0))
                      },
                      norm={
                          return(as.matrix(object@activities$normAct1))
                      },
                      threshold={
                          return(as.matrix(object@activities$thresholdAct2))
                      },
                      scaled={
                          return(as.matrix(object@activities$scaledAct3))
                      })
           })
})
