#' @rdname clinPredictionDenovo-methods
#' @aliases clinPredictionDenovo
setMethod("clinPredictionDenovo",
          signature=c(object="SigQuant"),
          definition=function(object, sampTrain, sigsTrain){
              message("this function is depreciated")
              # if(getExperiment(object)@feature.method != "drews"){
              #     stop("This function is only applicable to objects using drews method.")
              # }
              # # Load normalised signature activities
              # mNorm = object@activities$normAct1
              #
              # # Extract samples for training
              # if(is.null(sampTrain)) { stop("No sample names supplied.")}
              # if(is.null(sigsTrain)) { stop("No signature names supplied.")}
              # if(length(sigsTrain) != 2) { stop("So far only two signatures can be used.")}
              # mTrain = mNorm[ rownames(mNorm) %in% sampTrain, colnames(mNorm) %in% sigsTrain]
              # mTest = mNorm[ ! rownames(mNorm) %in% sampTrain, colnames(mNorm) %in% sigsTrain]
              #
              # # Scale training data and apply to test cohort
              # scaledTrain = scale(mTest)
              # lModel = list(mean = attr(scaledTrain, "scaled:center"),
              #               scale = attr(scaledTrain, "scaled:scale"))
              #
              # scaledTest = scaleByModel(mTest, lModel)
              #
              # # Do classification
              # vPred = ifelse(scaledTest[,sigsTrain[1]] >= scaledTest[,sigsTrain[2]], paste("Signature", sigsTrain[1], "higher"),
              #                paste("Signature", sigsTrain[2], "higher"))
              # return(vPred)
})
