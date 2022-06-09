#' getSamples
#'
#' Extracts sample names from a CNQuant object.
#'
#' @param object CNQuant object
#' @return A character vector
#' @export
#' @docType methods
#' @rdname getSamples-methods
#'
setGeneric("getSamples", function(object) standardGeneric("getSamples"))

#' getSegments
#'
#' Extracts copy number segment data from a CNQuant object.
#'
#' @param object CNQuant object
#' @return A data.frame
#' @export
#' @docType methods
#' @rdname getSegments-methods
#'
setGeneric("getSegments", function(object) standardGeneric("getSegments"))

#' getSamplefeatures
#'
#' Extracts sample feature data from a CNQuant object.
#'
#' @param object CNQuant object
#' @return A data.frame
#' @export
#' @docType methods
#' @rdname getSamplefeatures-methods
#'
setGeneric("getSamplefeatures", function(object) standardGeneric("getSamplefeatures"))

#' getSampleByComponent
#'
#' @param object CNQuant object
#' @return matrix containing the sample-by-component data
#' @export
#' @docType methods
#' @rdname getSampleByComponent-methods
#'
setGeneric("getSampleByComponent", function(object) standardGeneric("getSampleByComponent"))

#' getExperiment
#'
#' Extracts and returns copy number features from copy number profiles in a CNQuant object.
#'
#' @param object CNQuant object
#' @return A ExpQuant class object
#' @export
#' @docType methods
#' @rdname getExperiment-methods
#'
setGeneric("getExperiment", function(object) standardGeneric("getExperiment"))

#' calculateFeatures
#'
#' Extracts and returns copy number features from copy number profiles in a CNQuant object.
#'
#' @param object CNQuant object
#' @param method Method to extract copy number features. Default is "drews".
#' @param smooth.diploid Binary variable indicating whether segments close to 2 should be collapsed to 2 and merged together. Default is TRUE.
#' @param cores Number of CPU threads/cores to utilise via doParallel. Default is 1.
#' @return A CNQuant class object with extracted features stored in the "featData" slot
#' @export
#' @docType methods
#' @rdname calculateFeatures-methods
#'
setGeneric("calculateFeatures",function(object, method="drews", smooth.diploid=TRUE,cores=1)
    standardGeneric("calculateFeatures"))

#' calculateSampleByComponentMatrix
#'
#' Calculates and returns a sample-by-component matrix from copy number features in a CNQuant object.
#'
#' @param object CNQuant object
#' @param method Determines the mixture components used to calculate sum-of-posterior probabilities. Default is "drews".
#' @return A CNQuant class object with sum-of-posterior probabilities stored in the "featFitting" slot
#' @export
#' @docType methods
#' @rdname calculateSampleByComponentMatrix-methods
#'

setGeneric("calculateSampleByComponentMatrix",function(object, method="drews")
    standardGeneric("calculateSampleByComponentMatrix"))

#' calculateActivity
#'
#' Calculates and returns signature activities in a SigQuant object. Works best after function calculateSampleByComponentMatrix call. \cr \cr
#' The output of this function is a list of four matrices, the raw signature activities, the normalised activities, the normalised and
#' thresholded signature activities and the normalised, thresholded and scaled activities with the scaling factors obtained from the TCGA cohort.
#'
#'
#' @param object CNQuant object
#' @param method Determines the mixture components used to calculate sum-of-posterior probabilities. Default is "drews".
#' @return A SigQuant class object with four activity matrices stored in the "activities" slot
#' @export
#' @docType methods
#' @rdname calculateActivity-methods
#'

setGeneric("calculateActivity",function(object, method="drews")
    standardGeneric("calculateActivity"))

#' quantifyCNSignatures
#'
#' This function takes a copy number profile as input and returns signature activities.
#'
#' @param object CNQuant object
#' @param experimentName A user-specified name of the experiment
#' @param method The method used for calculating the signature activities. Default is "drews"
#' @param cores Number of threads/cores to use for parallel processing
#' @return A SigQuant class object with four activity matrices stored in the "activities" slot
#' @export
#' @docType methods
#' @rdname quantifyCNSignatures-methods
#'

setGeneric("quantifyCNSignatures",function(object, experimentName="Default", method="drews",cores=1)
    standardGeneric("quantifyCNSignatures"))

#' clinPredictionPlatinum
#'
#' The function takes signature activities based on Drews et al. methodology and predicts patient's response to platinum-based chemotherapies.
#'
#' @param object SigQuant object
#' @return A vector with "Predicted sensitive" or "Predicted resistant" for all samples in the input object.
#' @export
#' @docType methods
#' @rdname clinPredictionPlatinum-methods
#'

setGeneric("clinPredictionPlatinum",function(object)
    standardGeneric("clinPredictionPlatinum"))

#' clinPredictionDenovo
#'
#' The function takes signature activities based on Drews et al. methodology and predicts patient's response based on
#' user-specified pair of signatures. \cr \cr
#' The user should supply a vector of samples for training purposes. The function then trains the classifier on these samples before applying it to all samples and return the labels.
#'
#' @param object SigQuant object
#' @param sampTrain Vector of sample names that should be used for training the classifier.
#' @param sigsTrain Vector with two signature names on which the prediction should be based upon.
#' @return A vector with "Signature <1> higher" or "Signature <2> higher" for all samples in the input object.
#' @export
#' @docType methods
#' @rdname clinPredictionDenovo-methods
#'

setGeneric("clinPredictionDenovo",function(object, sampTrain, sigsTrain)
    standardGeneric("clinPredictionDenovo"))
