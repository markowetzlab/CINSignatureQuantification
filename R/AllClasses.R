#' ExpQuant class
#'
#' The `ExpQuant` class is an accessory class providing meta data documentation
#' used by `CNQuant` and `SigQuant` classes. The `ExpQuant` class stores:
#' Experiment name, initialisation time,last modified time, total sample count,
#' current sample count, genomic build, and specified method
#'
#' @slot experimentName Object of type `"character"` containing the name
#'   specified for the `CNQuant` class object by [createCNQuant()] or
#'   [quantifyCNSignatures()].
#' @slot init.date Object of type `"character"` containing the date at which
#'   copy number data was initially loaded into the `CNQuant` class object.
#' @slot last.modified Object of type `"character"` containing the date at which
#'   the `CNQuant` or `SigQuant` class object was last modified by a function
#'   which adds or modifies a slot of the described class.
#' @slot samples.full Object of type `"numeric"` containing total number of
#'   samples provided on object initialisation.
#' @slot samples.current Object of type `"numeric"` containing current number of
#'   samples within the object.
#' @slot build Object of type `"character"` containing the specified genome
#'   build (hg19 or hg38) by [createCNQuant()] or [quantifyCNSignatures()].
#' @slot feature.method Object of type `"character"` containing the specified
#'   copy number signatures method for the `CNQuant` class object by
#'   [createCNQuant()] or [quantifyCNSignatures()].
#' @importFrom methods new
#' @seealso [CNQuant-class]
#' @seealso [createCNQuant()]
#' @export
ExpQuant <- setClass("ExpQuant",
                     slots = list(experimentName = "character",
                                  init.date = "character",
                                  last.modified = "character",
                                  samples.full = "numeric",
                                  samples.current = "numeric",
                                  build = "character",
                                  feature.method = "character"
                     ),
                     prototype = list(experimentName = NULL,
                                      init.date = as.character(Sys.time()),
                                      last.modified = "NA",
                                      samples.full = NULL,
                                      samples.current = NULL,
                                      build="hg19",
                                      feature.method = "NA")
)

#' CNQuant class
#'
#' The `CNQuant` class is a structured S4 class object designed to contain copy
#' number and copy number signature-related data. This is initial class prior to
#' computation of copy number signatures which extends the `CNquant` class to
#' `SigQuant` class.
#'
#'
#' @slot segments Object of type `"list"` containing data.frame objects for each
#'   sample segment table. Access with [getSegments()].
#' @slot featData Object of type `"list"` containing data.frame objects for each
#'   copy number feature extracted by [calculateFeatures()]. Access with
#'   [getFeatures()].
#' @slot featFitting Object of type `"list"` containing objects feature
#'   component fitting, including the sample by component matrix and model
#'   parameters, calculated by [calculateSampleByComponentMatrix()]. Access with
#'   [getSampleByComponent()].
#' @slot samplefeatData Object of type `"data.frame"` containing sample-level
#'   information related to samples provided on `CNQuant` class initialisation,
#'   including segment counts and sample ploidy. Access with
#'   [getSamplefeatures()] and add additional sample features using [addsampleFeatures]
#' @slot ExpData An object of class `"ExpQuant"`, see [ExpQuant-class].
#' @importFrom methods new
#' @seealso [SigQuant-class]
#' @seealso [createCNQuant()]
#' @export
CNQuant <- setClass("CNQuant",
                    slots = list(segments = "list",
                                 featData = "list",
                                 featFitting = "list",
                                 samplefeatData = "data.frame",
                                 ExpData = "ExpQuant")
)

#' SigQuant class
#'
#' The `SigQuant` class is an extension of the `CNQuant` S4 class object
#' designed to contain copy number and copy number signature-related data. This
#' class is initialised on generation of copy number signature definitions or
#' activities from `CNQuant` class object with the necessary fitted feature
#' data. It inherits all slots from `CNQuant` and adds six additional slots.
#'
#' @slot activities Object of type `"list"` containing matrices for each
#'   calculated copy number signature activity matrix including raw, normalised,
#'   threshold-adjusted, and scaled signature activities. Access with
#'   [getActivities()].
#' @slot signature.model Object of type `"character"` specifying the methodology
#'   used to calculate copy number signatures.
#' @slot backup.signatures bject of type `"matrix"` containing the copy number
#'   signature defintions matrix (signature-by-component) used to calculate the
#'   copy number signature activities.
#' @slot backup.thresholds Object of type `"numeric"` containing the uniform or
#'   signature-specific thresholds applied to the normalised copy number
#'   signature activities at which to reduce activity to zero based on signature
#'   noise simulations.
#' @slot backup.scale Object of type `"list"` containing two numeric vectors
#'   corresponding to the mean and scaling parameter for each of the copy number
#'   signature activities used to generate the scaled activities matrix.
#' @slot backup.scale.model Object of type `"character"` specifying the scale
#'   model used to generate the scaled signature activities.
#' @importFrom methods new
#' @seealso [CNQuant-class]
#' @seealso [createCNQuant()]
#' @export
SigQuant <- setClass("SigQuant",
                              contains = "CNQuant",
                              slots = list(
                                  activities = "list",
                                  signature.model = "character",
                                  backup.signatures = "matrix",
                                  backup.thresholds = "numeric",
                                  backup.scale = "list",
                                  backup.scale.model = "character"
                              )
)
