#' @importFrom limSolve lsei
LinCombDecompSigs = function (component_by_sample, component_by_signature) {

    # Failsafe: Convert signatures to matrix
    mSignatures = as.matrix(component_by_signature)

    # Prepare LCD (needs matrix and vectors for each signature)
    numSigs = ncol(mSignatures)
    G = diag(numSigs)
    H = rep(0, numSigs)

    # Prepare output
    dfOutExp = data.frame()

    # Loop over signatures
    numSamps = ncol(component_by_sample)
    for (i in seq_len(numSamps)) {

        # Perform the magic. LCD on a signature (vector of weights) and the input matrix
        lLCDResult = lsei(A = mSignatures,
                                    B = component_by_sample[, i],
                                    G = G,
                                    H = H,
                                    verbose = FALSE)

        # Extract relevant vector
        vExp = as.vector(lLCDResult$X)

        # Attach results to output data frame
        dfOutExp[seq(1, numSigs, 1), i] = vExp

        # Clean up
        rm(lLCDResult)
    }

    # Transfer names and return
    colnames(dfOutExp) = colnames(component_by_sample)
    rownames(dfOutExp) = colnames(component_by_signature)

    return(dfOutExp)
}
