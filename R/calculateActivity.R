#' @rdname calculateActivity-methods
#' @aliases calculateActivity
setMethod("calculateActivity",
          signature=c(object="CNQuant"),
          definition=function(object, method=NULL,cancer.subset=NULL){
              if(length(object@featFitting) == 0){
                  stop("Sample-by-component unavailable - run 'calculateSampleByComponentMatrix()'")
              }
              # Check method
              if(is.null(method)){
                  method <- getExperiment(object)@feature.method
              }
              switch(method,
                     mac={
                         if(!is.null(cancer.subset)){
                             stop("argument 'cancer.subset = '",cancer.subset,"' provided - cancer subsetting not available for mac method",call. = F)
                         }
                         SigActs <- calculateActivityMac(object)
                         Hraw <- t(SigActs[[1]])
                         SigActs <- t(SigActs[[2]])
                         SigActs <- list(rawAct0=Hraw,normAct1=NULL,thresholdAct2=SigActs,scaledAct3=NULL)

                         #W<-t(get(load("data/Macintyre2018_OV_Signatures_normalised.rda")))
                         W <- t(get(utils::data("Macintyre2018_OV_Signatures_normalised",envir = environment())))
                         # Combine results
                         methods::new("SigQuant",object,
                                      activities=SigActs,
                                      signature.model = method,
                                      backup.signatures=W,
                                      backup.thresholds=0.01,
                                      backup.scale=list(mean=c("NULL"),weight=c("NULL")),
                                      backup.scale.model="NULL",
                                      ExpData = methods::initialize(object@ExpData,
                                                                    last.modified = as.character(Sys.time()),
                                                                    feature.method = method))
                     },
                     drews={
                         # Calculate activities
                         Hraw = calculateActivityDrews(object,cancer.subset=cancer.subset)
                         # Apply thresholds, normalisation and TCGA-specific scaling
                         lSigs = applyThreshNormAndScaling(Hraw,cancer.subset=cancer.subset)

                         # Load data to be put into model as backup
                         #W = get(load("data/Drews2022_TCGA_Signatures.rda"))
                         W = get(utils::data("Drews2022_TCGA_Signatures",envir = environment()))
                         #vThresh = get(load("data/Drews2022_TCGA_Signature_Thresholds.rda"))
                         vThresh = get(utils::data("Drews2022_TCGA_Signature_Thresholds",envir = environment()))
                         #lScales = get(load("data/Drews2022_TCGA_Scaling_Variables.rda"))
                         lScales = get(utils::data("Drews2022_TCGA_Scaling_Variables",envir = environment()))
                         # Combine results
                         methods::new("SigQuant",object,
                                      activities=lSigs,
                                      signature.model = method,
                                      backup.signatures=W,
                                      backup.thresholds=vThresh,
                                      backup.scale=lScales,
                                      backup.scale.model="TCGA",
                                      ExpData = methods::initialize(object@ExpData,
                                                                    last.modified = as.character(Sys.time()),
                                                                    feature.method = method))

                     })
              })
