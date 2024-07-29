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
              } else if(method != getExperiment(object)@feature.method){
                  featMethod <- getExperiment(object)@feature.method
                  sxcMethodStop <- paste0("provided method different to feature method - provided: ",method," | feature method: ",featMethod)
                  stop(sxcMethodStop)
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
