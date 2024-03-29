startCopynumberFeatureExtractionDrews = function(dtSmooth, cores = 1, RMNORM = TRUE, build="hg19") {

    # Convert to data frame
    dfBR = data.frame(dtSmooth)
    # Split by sample
    lBR = split( dfBR, dfBR$sample )
    # Extract features
    brECNF = extractCopynumberFeaturesDrews(lBR, cores = cores, rmNorm = RMNORM, build=build)

    return(brECNF)
}
