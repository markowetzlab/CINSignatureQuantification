## Function for identifying and merging neighbouring segments that are close to each other (as defined by a user-supplied threshold)
smoothAndMergeSegments = function(dfAllSegs, CORES, WIGGLE = 0.1, colNameSegVal = "segVal", colNameChr = "chromosome", IGNOREDELS = FALSE) {

    # Explicit conversion to numeric. Just to be on the safe site.
    dfAllSegs$start = as.numeric( dfAllSegs$start )
    dfAllSegs$end = as.numeric( dfAllSegs$end )
    dfAllSegs$segVal = as.numeric( dfAllSegs$segVal )

    # Set everything very close to 2 to 2
    dfAllSegs$segVal[ dfAllSegs$segVal > (2-WIGGLE) & dfAllSegs$segVal < (2+WIGGLE) ] = 2
    # Merge segments only when two normal follow each other -> SMOOTHINGFACTOR = 0
    dfAllSegs = idSmoothingTargets(dfAllSegs, WIGGLE = 0, colNameSegVal = "segVal", colNameChr = "chromosome", IGNOREDELS = IGNOREDELS)
    # Split by sample name
    lRaw = split(dfAllSegs, dfAllSegs$sample)

    # Smooth segments by taking the weighted average of the segVal and their lengths
    lSmooth = smoothSegments(lRaw, CORES, SMOOTHINGFACTOR = 0, colNameMerge = "segVal", colNameChr = "chromosome",
                             colNameStart = "start", colNameEnd = "end", IGNOREDELS = IGNOREDELS, asDf = FALSE)
    dtSmooth = rbindlist(lSmooth)

    return(dtSmooth)
}
