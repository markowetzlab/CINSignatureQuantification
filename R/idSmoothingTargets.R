# Function for identifying neighbouring segments that have the same copy number and should be merged into one segment
idSmoothingTargets = function(dfAllSegs, WIGGLE, colNameSegVal, colNameChr, IGNOREDELS = TRUE) {

    ### Check column name
    testSegVal = dfAllSegs[[colNameSegVal]][1]
    testChr = dfAllSegs[[colNameChr]][1]

    # Quick sanity checks
    if(! is.numeric(testSegVal)) { stop("Segment Value column has no numeric value in it. Supplied correct column name? Forgot conversion?")}
    if(is.null(testSegVal)) { stop("Chromosome column has no numeric value in it. Supplied correct column name?")}

    # Take differences to segment down below
    dfAllSegs$diffs = c( abs( dfAllSegs[[colNameSegVal]][1:(nrow(dfAllSegs)-1)] - dfAllSegs[[colNameSegVal]][2:nrow(dfAllSegs)] ), WIGGLE+1)
    # Set TRUE if difference to next segment is smaller than the user supplied cutoff
    dfAllSegs$smooth = dfAllSegs$diffs <= WIGGLE
    # Set all segments which are last in a chromosome to FALSE. This also prevents leaking to other samples and cohorts.
    dfAllSegs$smooth[ cumsum( rle(as.character(dfAllSegs[[colNameChr]]))$lengths ) ] = FALSE

    # Ignore deletions if wished
    if(IGNOREDELS) { dfAllSegs$smooth[ dfAllSegs[[colNameSegVal]] == 0 ] = FALSE }

    return( dfAllSegs )
}
