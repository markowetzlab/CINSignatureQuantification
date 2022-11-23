applyThreshNormAndScaling = function(Hraw,cancer.subset=NULL) {

    # Normalise matrix
    H = t( apply(Hraw, 2, function(x) x/sum(x)) )

    # Apply signature-specific thresholds (no renormalising happening in order to avoid inflation of signal)
    #vThresh = get(load("data/Drews2022_TCGA_Signature_Thresholds.rda"))
    vThresh = get(utils::data("Drews2022_TCGA_Signature_Thresholds",envir = environment()))
    if(!is.null(cancer.subset)){
        subset = getCancerSpecificSignatures(cancer.subset)
        vThresh = vThresh[names(vThresh) %in% subset]
    }
    threshH = sapply(names(vThresh), function(thisSig) {

        sigVals = H[,thisSig]
        sigVals[ sigVals < vThresh[thisSig] ] = 0

        return(sigVals)
    })

    if(is.null(dim(threshH))){
        threshH <- t(as.matrix(threshH))
        rownames(threshH) <- rownames(H)
    }
    # Scale according to TCGA-specific scaling factors
    #lScales = get(load("data/Drews2022_TCGA_Scaling_Variables.rda"))
    lScales = get(utils::data("Drews2022_TCGA_Scaling_Variables",envir = environment()))
    if(!is.null(cancer.subset)){
        subset = getCancerSpecificSignatures(cancer.subset)
        lScales = lapply(lScales,FUN=function(x){x[names(x) %in% subset]})
    }

    threshScaledH = scaleByModel(threshH, lScales)

    # Combine for return
    lOut = list(rawAct0 = t(Hraw), normAct1 = H, thresholdAct2 = threshH, scaledAct3 = threshScaledH)
    return(lOut)
}
