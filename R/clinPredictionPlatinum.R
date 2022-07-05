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
              lModel = get(utils::data("Drews2022_CX3CX2_Clinical_classifier",envir = environment()))
              if(is.null(dim(mNorm[,names(lModel$mean)]))){
                  H <- rownames(mNorm)
                  mNorm <- mNorm[,names(lModel$mean)]
                  mNorm <- t(as.matrix(mNorm))
                  rownames(mNorm) <- H
                  mNormGBRCA1 = scaleByModel(mNorm, lModel)
              } else {
                  mNormGBRCA1 = scaleByModel(mNorm[,names(lModel$mean)], lModel)
              }
              # Do classification
              vPred = ifelse(mNormGBRCA1[,"CX3"] >= mNormGBRCA1[,"CX2"], "Predicted sensitive", "Predicted resistant")
              return(vPred)
})
