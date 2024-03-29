#' @rdname calculateFeatures-methods
#' @aliases calculateFeatures
#' @importFrom data.table rbindlist
setMethod("calculateFeatures",
          signature=c(object="CNQuant"),
          definition=function(object,method = NULL,smooth.diploid=TRUE,cores=1,DCIN = 20){
              methods <- c("mac","drews")
              if(is.null(method)){
                  stop("no method provided")
              } else if(!method %in% methods){
                  stop("unknown method")
              }
              switch(method,
                mac={
                    featData <- extractCopynumberFeaturesMac(object@segments,cores = cores,
                    build=object@ExpData@build)
                    methods::initialize(object,featData=featData,
                           ExpData = methods::initialize(object@ExpData,
                                                         last.modified = as.character(Sys.time()),
                                                         feature.method = method))
                },
                drews={
                    # Smooth and merge segments
                    smoothed = smoothAndMergeSegments(getSegments(object),
                                                      CORES = cores,
                                                      WIGGLE = 0.1,
                                                      colNameSegVal = "segVal",
                                                      colNameChr = "chromosome",
                                                      IGNOREDELS = FALSE)
                    # Avoid measurement errors
                    smoothed <- avoidMeasurementErrors(smoothed)
                    # Filter samples without CIN
                    filtered <- removeQuietSamples(smoothed, DCIN = DCIN)

                    kept <- ifelse(unique(smoothed$sample) %in% unique(filtered$sample),TRUE,FALSE)
                    if(any(!kept)){
                        message(cat("Low segment counts: Features not computed for",
                                    length(unique(smoothed$sample)[!kept]),"of",
                                    length(unique(smoothed$sample)),
                                           "samples - see `getSampleFeatures()`"))
                    }
                    object <- addsampleFeatures(object = object,
                                                sample.data = data.frame(sample=unique(smoothed$sample),
                                                                         computeSigs=kept))
                    if(nrow(filtered) == 0){
                        warning("no samples with sufficient copy number alteration counts")
                        methods::initialize(object,featData=list(),
                                            ExpData = methods::initialize(object@ExpData,
                                                                          last.modified = as.character(Sys.time()),
                                                                          feature.method = method))
                    } else {
                    # Extract
                        featData = startCopynumberFeatureExtractionDrews(filtered, cores = cores, RMNORM = TRUE, build=object@ExpData@build)
                        methods::initialize(object,featData=featData,
                                        ExpData = methods::initialize(object@ExpData,
                                                                      last.modified = as.character(Sys.time()),
                                                                      feature.method = method))
                    }
                })
          })
