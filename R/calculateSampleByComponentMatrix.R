#' @rdname calculateSampleByComponentMatrix-methods
#' @aliases calculateSampleByComponentMatrix
setMethod("calculateSampleByComponentMatrix",
          signature=c(object="CNQuant"),
          definition=function(object, method=NULL){
              if(length(object@featData) == 0){
                  stop("Copy number features unavailable (run 'calculateFeatures()) or could not be calculated due to too few segments'")
              }
              if(is.null(method)){
                method <- getExperiment(object)@feature.method
              }
              switch(method,
                     mac={
                         sxc <- calculateSampleByComponentMatrixMac(object@featData,
                                                                      UNINFPRIOR = FALSE)
                         sxc = c(sxc, method = method)
                         methods::new("CNQuant",object,featFitting=sxc,
                                      ExpData = methods::initialize(object@ExpData,
                                                                    last.modified = as.character(Sys.time()),
                                                                    feature.method = method))
                     },
                     drews={
                         lSxC = calculateSampleByComponentMatrixDrews(object@featData,
                                                                     UNINFPRIOR = TRUE)
                         lSxC = c(lSxC, method = method)
                         methods::new("CNQuant",object,featFitting=lSxC,
                                      ExpData = methods::initialize(object@ExpData,
                                                                    last.modified = as.character(Sys.time()),
                                                                    feature.method = method))
                     })

          })
