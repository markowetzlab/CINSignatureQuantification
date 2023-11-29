calculateActivityDrews = function(object,cancer.subset=NULL) {

    # Extract relevant information from object
    V = object@featFitting$sampleByComponent
    nSamp = nrow(object@featFitting$sampleByComponent)
    nFeat = ncol(object@featFitting$sampleByComponent)

    # Load signatures
    #W = get(load("data/Drews2022_TCGA_Signatures.rda"))
    W = get(utils::data("Drews2022_TCGA_Signatures",envir = environment()))
    if(!is.null(cancer.subset)){
        subset = getCancerSpecificSignatures(cancer.subset)
        W = W[rownames(W) %in% subset,]
    }
    # Sanity check mutational catalogue (not really necessary)
    if(nSamp > nFeat) {
        # Case 1: More samples than features
        if(nrow(V) > ncol(V)) { V = t(V) }
    } else if(nSamp < nFeat) {
        # Case 2: Fewer samples than features
        if(nrow(V) < ncol(V)) { V = t(V) }
    } else {
        # Case 3: Edge case where there are as many samples as features
        if(sum(grepl("segsize", colnames(V))) > 0) { V = t(V) }
    }


    # Sanity check signature matrix
    if(nrow(W) < ncol(W)) W = t(W)
    # Check order of components and fix if necessary
    if(! identical(rownames(W), rownames(V)) ) {
        W = W[ match(rownames(V), rownames(W)), ]
    }

    ### Functions needs:
    ## Full matrix V        mutCatalogue        components (rows) by samples (cols)    <= HAVE
    ## Left matrix W        sigCatalogue        components (rows) by signature (cols)   <= HAVE
    ## Right matrix H       expCatalogue        signature (rows) by samples (cols)     <= WANT

    # component_by_sample => NxM => N - Features, M - Samples => Component by Sample matrix
    # component_by_signature => NxL => N - Features, L - Signatures => Component by Signature matrix
    Hraw = as.matrix(LinCombDecompSigs(component_by_sample = V, component_by_signature = W))
    return(Hraw)

}
