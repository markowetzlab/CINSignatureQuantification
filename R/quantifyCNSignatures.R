#' @rdname quantifyCNSignatures-methods
#' @aliases quantifyCNSignatures
setMethod(
    "quantifyCNSignatures",
    signature = c(object = "ANY"),
    definition = function(object,
                          experimentName = "Default",
                          method = "drews",
                          cores = 1,
                          build = "hg19",
                          cancer.subset = NULL) {
        # Check method
        if (is.null(method) | !(method %in% c("drews", "mac"))) {
            stop("Method was neither 'drews' nor 'mac'.")
        }

        # Create object from CN profiles
        # TODO: Extend for QDNAseq
        cigTCGA = createCNQuant(data = object, experimentName = experimentName,
                                build = build)
        # Extract features
        cigTCGA = calculateFeatures(object = cigTCGA,
                                    method = method,
                                    cores = cores)
        # Calculate sum-of-posterior matrix
        cigTCGA = calculateSampleByComponentMatrix(object = cigTCGA, method =
                                                       method)
        # Calculate signature activities
        cigTCGA = calculateActivity(object = cigTCGA, method = method,
                                    cancer.subset = cancer.subset)
        return(cigTCGA)
    }
)
