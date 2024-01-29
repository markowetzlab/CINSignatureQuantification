#' getSamples
#'
#' Returns sample names from a CNQuant class object.
#'
#' @param object CNQuant object
#' @return character vector containing all available sample names
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- createCNQuant(TCGA_478_Samples_SNP6_GOLD)
#'   Samples <- getSamples(cnobj[1:10])
#' @export
#' @docType methods
#' @rdname getSamples-methods
#'
setGeneric("getSamples", function(object) standardGeneric("getSamples"))

#' getSampleByComponent
#'
#' Returns the sample-by-component a `CNQuant` or `SigQuant` class object.
#' @param object CNQuant object with fitted feature data
#' @return matrix containing the sample-by-component data
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- createCNQuant(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- calculateFeatures(cnobj[1:10],method="drews")
#'   cnobj <- calculateSampleByComponentMatrix(cnobj)
#'   SampleByComponent <- getSampleByComponent(cnobj)
#' @export
#' @docType methods
#' @rdname getSampleByComponent-methods
#'
setGeneric("getSampleByComponent", function(object) standardGeneric("getSampleByComponent"))

#' getSamplefeatures
#'
#' Extracts sample feature data from a `CNQuant` class object.
#'
#' @param object CNQuant object
#' @return A data.frame containing the sample-specific features (ploidy, segment
#'   count, etc.)
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- createCNQuant(TCGA_478_Samples_SNP6_GOLD)
#'   SampleFeats <- getSamplefeatures(cnobj[1:10])
#' @seealso [addsampleFeatures()]
#' @export
#' @docType methods
#' @rdname getSamplefeatures-methods
#'
setGeneric("getSamplefeatures", function(object) standardGeneric("getSamplefeatures"))

#' getSegments
#'
#' Extracts copy number segment data from a `CNQuant` class object.
#'
#' @param object CNQuant object
#' @return A data.frame containing the segment table for all samples.
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- createCNQuant(TCGA_478_Samples_SNP6_GOLD)
#'   segments <- getSegments(cnobj[1:10])
#' @export
#' @docType methods
#' @rdname getSegments-methods
#'
setGeneric("getSegments", function(object) standardGeneric("getSegments"))

#' getActivities
#'
#' Extracts signature activities from a `SigQuant` class object.
#'
#' @param object SigQuant object
#' @param type type of copy number signature matrix to return 'raw', 'norm',
#'   'threshold', or 'scaled' (default: threshold).
#' @details Type parameter can be specified to return a specific type of
#'   activities matrix. For most users, the threshold-corrected ("threshold")
#'   activities would be preferred and are returned by default if no type is
#'   specified. * `raw`: Copy number signature activities without any
#'   normalisation or zero-exposure correction. * `normalised`: Copy number
#'   signature activities normalised to between 0 and 1.0 where each sample sums
#'   to 1.0. * `threshold`: Copy number signature activities which are
#'   normalised and have signature activities less than a given threshold
#'   reduced to zero. This threshold application is different for each method
#'   where `drews` applies a signature-specific threshold and `mac` applies a
#'   uniform threshold for all signatures. Threshold-corrected signatures using
#'   `drews` method may not sum to 1.0 as normalisation is not reapplied so as
#'   to not inflate signatures artificially. * `scaled`: Copy number signature
#'   activities scaled by the signature activity distributions identified in
#'   drews et al. 2022 which aims to make copy number signature activities from
#'   different data sets comparable to those generated in the drews et al. 2022
#'   publication.
#' @return A numeric matrix containing copy number signature activities
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- createCNQuant(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- calculateFeatures(cnobj[1:10],method="drews")
#'   cnobj <- calculateSampleByComponentMatrix(cnobj)
#'   cnobj <- calculateActivity(cnobj)
#'   activities <- getActivities(cnobj,type="threshold")
#' @export
#' @docType methods
#' @rdname getActivities-methods
#'
setGeneric("getActivities", function(object,type="threshold") standardGeneric("getActivities"))

#' getExperiment
#'
#' Returns the experiment details from a `CNQuant` class object.
#'
#' @param object CNQuant object
#' @return A `ExpQuant` class object
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- createCNQuant(TCGA_478_Samples_SNP6_GOLD)
#'   getExperiment(cnobj)
#' @seealso [ExpQuant-class]
#' @export
#' @docType methods
#' @rdname getExperiment-methods
#'
setGeneric("getExperiment", function(object) standardGeneric("getExperiment"))

#' getFeatures
#'
#' Returns the copy number features from a `CNQuant` class object. Features
#' available will be those extracted as defined by the copy number signature
#' method used.
#'
#' @param object CNQuant object containing copy number features
#' @param feat name or vector names of features to return. By default returns
#'   all features.
#' @return list containing the extracted feature data
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- createCNQuant(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- calculateFeatures(cnobj[1:10],method="drews")
#'   features <- getFeatures(cnobj)
#' @export
#' @docType methods
#' @rdname getFeatures-methods
#'
setGeneric("getFeatures", function(object,feat=NULL) standardGeneric("getFeatures"))

#' calculateFeatures
#'
#' Calculates and returns copy number features from copy number profiles in a
#' `CNQuant` class object. The output from this function is stored within the
#' `featData` slot in the returned `CNQuant` class object.
#'
#' @param object CNQuant object
#' @param method Method to extract copy number features. Default is "drews".
#' @param smooth.diploid Binary variable indicating whether segments close to 2
#'   should be collapsed to 2 and merged together. Default is TRUE.
#' @param cores Number of CPU threads/cores to utilise via doParallel. Default
#'   is 1. Maximum number is equal to the number of features to extract (drews &
#'   mac methods = 6 features).
#' @param DCIN Threshold for required number of non-diploid segments to compute
#'   copy number features (and subsequently copy number signatures) using method
#'   "drews". Default is 20. This parameter should not need to be changed and
#'   will affect feature values and signature activity.
#' @return A CNQuant class object with extracted features stored in the
#'   "featData" slot
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- createCNQuant(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- calculateFeatures(cnobj[1:10],method="drews")
#' @seealso [getFeatures()]
#' @export
#' @docType methods
#' @rdname calculateFeatures-methods
#'
setGeneric("calculateFeatures",function(object, method="drews",smooth.diploid=TRUE,cores=1,DCIN = 20)
    standardGeneric("calculateFeatures"))

#' calculateSampleByComponentMatrix
#'
#' Calculates and returns a sample-by-component matrix from copy number features
#' in a `CNQuant` class object. The output from this function is stored within
#' the `featFitting` slot in the returned `CNQuant` class object.
#'
#' @param object CNQuant object
#' @param method Determines the mixture components used to calculate
#'   sum-of-posterior probabilities. Default is "drews".
#' @return A CNQuant class object with sum-of-posterior probabilities stored in
#'   the "featFitting" slot
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- createCNQuant(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- calculateFeatures(cnobj[1:10],method="drews")
#'   cnobj <- calculateSampleByComponentMatrix(cnobj)
#' @seealso [getSampleByComponent()]
#' @export
#' @docType methods
#' @rdname calculateSampleByComponentMatrix-methods
#'
setGeneric("calculateSampleByComponentMatrix",function(object, method="drews")
    standardGeneric("calculateSampleByComponentMatrix"))

#' calculateActivity
#'
#' Calculates and returns signature activities in a `SigQuant` class object
#' from Sample-By-Component data.
#'
#' @param object CNQuant object
#' @param method Determines the mixture components used to calculate
#'   sum-of-posterior probabilities. Default is "drews".
#' @param cancer.subset A TCGA cancer subset (e.g. BRCA or KIRC) identifier to
#'   select a subset of signatures for calculating signature activity.
#' @return A SigQuant class object with four activity matrices stored in the
#'   "activities" slot
#' @details The output of this function is a list of four matrices, the raw
#'   signature activities, the normalised activities, the normalised and
#'   thresholded signature activities and the normalised, thresholded and scaled
#'   activities with the scaling factors obtained from the TCGA cohort.
#'   * `raw`:
#'   Copy number signature activities without any normalisation or zero-exposure
#'   correction.
#'   * `normalised`: Copy number signature activities normalised to
#'   between 0 and 1.0 where each sample sums to 1.0.
#'   * `threshold`: Copy number
#'   signature activities which are normalised and have signature activities
#'   less than a given threshold reduced to zero. This threshold application is
#'   different for each method where `drews` applies a signature-specific
#'   threshold and `mac` applies a uniform threshold for all signatures.
#'   Threshold-corrected signatures using `drews` method may not sum to 1.0 as
#'   normalisation is not reapplied so as to not inflate signatures
#'   artificially.
#'   * `scaled`: Copy number signature activities scaled by the
#'   signature activity distributions identified in drews et al. 2022 which aims
#'   to make copy number signature activities from different data sets
#'   comparable to those generated in the drews et al. 2022 publication.
#' @return A numeric matrix containing copy number signature activities
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- createCNQuant(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- calculateFeatures(cnobj[1:10],method="drews")
#'   cnobj <- calculateSampleByComponentMatrix(cnobj)
#'   cnobj <- calculateActivity(cnobj)
#' @seealso [getActivities()]
#' @export
#' @docType methods
#' @rdname calculateActivity-methods
#'
setGeneric("calculateActivity",function(object, method="drews", cancer.subset=NULL)
    standardGeneric("calculateActivity"))

#' quantifyCNSignatures
#'
#' This function acts a wrapper function for the entire copy number signature
#' pipeline without the need to specify each step. It takes a copy number
#' profile as input and returns a `SigQuant` class object containing calculated
#' copy number signature activities and all required preceding data.
#'
#' @param object input data
#' @param experimentName A user-specified name of the experiment
#' @param method The method used for calculating the signature activities.
#'   Default is "drews"
#' @param cores Number of threads/cores to use for parallel processing
#' @param build Genome build to use, either hg19 or hg38 (default: hg19)
#' @param cancer.subset A TCGA cancer subset (e.g. BRCA or KIRC) identifier to
#'   select a subset of signatures for calculating signature activity
#' @return A SigQuant class object with four activity matrices stored in the
#'   "activities" slot
#' @details
#'   * object: Input for the object parameter should be either a segment
#'   table, with the fields chromosome, start, end, segVal, and sample or a
#'   `QDNAseqCopyNumbers` class object from the [QDNAseq
#'   package](https://github.com/ccagc/QDNAseq). Segment table input can be an
#'   existing `data.frame` or a file path to an appropriately delimited file
#'   containing the segment table data.
#'   * experimentName: experimentName can be
#'   a character string to name the `CNQuant` class object for future reference.
#'   It currently has no usage in any functions.
#'   * method: method can be either
#'   "drews" or "mac" and is saved within the Experiment metadata. It is
#'   utilised in downstream functions to implement a consistent methodology.
#'   * cores: Multi-thread/multi-core processing is implemented for a number of
#'   functions to decrease computation time and is implemented through the
#'   [doParallel
#'   package](https://cran.r-project.org/web/packages/doParallel/index.html).
#'   * build: character string to specify the genome build to use when extracting
#'   copy number features. Only human data using either hg19 or hg38 is
#'   currently supported.
#'   * cancer.subset: cancer.subset parameter allows users
#'   to calculate copy number signature activities for subsets of the copy
#'   number signatures defined in drews et al. 2022 which are specific to a
#'   given cancer type. Cancer subsets match the initialism used by TCGA (e.g
#'   KIRC, BRCA).
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   t478 <- TCGA_478_Samples_SNP6_GOLD
#'   subsample <- unique(t478$sample)[1:10]
#'   t478 <- t478[t478$sample %in% subsample]
#'   cnobj <- quantifyCNSignatures(t478,method="drews")
#' @seealso [createCNQuant()]
#' @export
#' @docType methods
#' @rdname quantifyCNSignatures-methods
#'
setGeneric("quantifyCNSignatures",function(object, experimentName="Default",
                                           method="drews",cores=1,
                                           build="hg19",cancer.subset=NULL)
    standardGeneric("quantifyCNSignatures"))

#' clinPredictionPlatinum
#'
#' The function takes signature activities based on Drews et al. 2022
#' methodology and predicts response to platinum-based chemo-therapeutics.
#'
#' @param object SigQuant object
#' @return A vector with "Predicted sensitive" or "Predicted resistant" for all
#'   samples in the input object.
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- createCNQuant(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- calculateFeatures(cnobj[1:10],method="drews")
#'   cnobj <- calculateSampleByComponentMatrix(cnobj)
#'   cnobj <- calculateActivity(cnobj)
#'   cnobj <- clinPredictionPlatinum(cnobj)
#' @export
#' @docType methods
#' @rdname clinPredictionPlatinum-methods
#'
setGeneric("clinPredictionPlatinum",function(object)
    standardGeneric("clinPredictionPlatinum"))

#' clinPredictionDenovo
#'
#' The function takes signature activities based on Drews et al. methodology and
#' predicts patient's response based on user-specified pair of signatures. The
#' user should supply a vector of samples for training purposes. The function
#' then trains the classifier on these samples before applying it to all samples
#' and return the labels.
#'
#' @param object SigQuant object
#' @param sampTrain Vector of sample names that should be used for training the
#'   classifier.
#' @param sigsTrain Vector with two signature names on which the prediction
#'   should be based upon.
#' @return A vector with "Signature <1> higher" or "Signature <2> higher" for
#'   all samples in the input object.
#' @examples
#'   data(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- createCNQuant(TCGA_478_Samples_SNP6_GOLD)
#'   cnobj <- calculateFeatures(cnobj[1:10],method="drews")
#'   cnobj <- calculateSampleByComponentMatrix(cnobj)
#'   cnobj <- calculateActivity(cnobj)
#'   # Not a real predictor
#'   cnobj <- clinPredictionDenovo(cnobj,
#'              sampTrain = getSamples(cnobj[1:5]),
#'              sigsTrain = c("CX1","CX6"))
#' @export
#' @docType methods
#' @rdname clinPredictionDenovo-methods
#'
setGeneric("clinPredictionDenovo",function(object, sampTrain, sigsTrain)
    standardGeneric("clinPredictionDenovo"))


#' exportCIN
#'
#' This function exports data from a `CNQuant` or `SigQuant` class object.
#'
#' @param object CNQuant or SigQuant object
#' @param outputDir Output location provided as a writable file directory. This
#'   defaults to the current working directory.
#' @param outputPrefix A prefix added to the beginning of exported files.
#' @param sep Default file seperator used in writing files (Default: '\t')
#' @param fullExport Provide a full export of all data contained within the
#'   provided object, including reference data, feature values, models, and
#'   supporting information (default: FALSE).
#'
#' @return NULL
#' @examples
#' @export
#'
#' @docType methods
#' @rdname exportCIN-methods
#'
setGeneric("exportCIN",function(object,outputDir=NULL,outputPrefix=NULL,sep="\t",fullExport=FALSE)
    standardGeneric("exportCIN"))
