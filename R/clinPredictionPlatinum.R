#' @rdname clinPredictionPlatinum-methods
#' @aliases clinPredictionPlatinum
setMethod("clinPredictionPlatinum",
          signature=c(object="SigQuant"),
          definition=function(object){
              if(getExperiment(object)@feature.method != "drews"){
                  stop("This function is only applicable to objects using drews method.")
              }
              # Load normalised signature activities
              mNorm = object@activities$normAct1

              # Load and apply gBRCA1 scaling vars
              lModel = get(load("data/Drews2022_CX3CX2_Clinical_classifier.rda"))
              mNormGBRCA1 = scaleByModel(mNorm[,names(lModel$mean)], lModel)

              # Do classification
              vPred = ifelse(mNormGBRCA1[,"CX3"] >= mNormGBRCA1[,"CX2"], "Predicted sensitive", "Predicted resistant")
              return(vPred)
})
