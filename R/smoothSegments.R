#' @importFrom data.table setDT
## Smooth segments that are close to
smoothSegments = function(lRaw, CORES, WIGGLE, colNameMerge, colNameChr, colNameStart, colNameEnd,
                          IGNOREDELS = TRUE, asDf = FALSE) {

    ### Check column names
    test = lRaw[[1]]
    testMerge = test[[colNameMerge]][1]
    testChr = test[[colNameChr]][1]
    testStart = test[[colNameStart]][1]
    testEnd = test[[colNameEnd]][1]
    if(! is.numeric(testMerge)) { stop("Merge column has no numeric value in it. Supplied correct column name?")}
    if(is.null(testChr)) { stop("Chromosome column has no numeric value in it. Supplied correct column name?")}
    if(! is.numeric(testStart)) { stop("Start column has no numeric value in it. Supplied correct column name?")}
    if(! is.numeric(testEnd)) { stop("End column has no numeric value in it. Supplied correct column name?")}

    # Add diff column to names we want to keep when merging (comes from function "idSmoothingTargets").
    colNameMerge = c(colNameMerge, "diffs")
    if(CORES > 1) {
        if (!requireNamespace("doParallel", quietly = TRUE)) {
            stop(
                "Package \"doParallel\" must be installed to use multiple threads/cores.",
                call. = FALSE
            )
        }
        `%dopar%` <- foreach::`%dopar%`
        doParallel::registerDoParallel(CORES)
        lSmooth = foreach::foreach(thisSample = lRaw, .final = function(x) stats::setNames(x, names(lRaw)) ) %dopar% {

            thisOut = thisSample
            stillSmoothing = sum(thisOut$smooth)
            while( stillSmoothing > 0 ) {
                # For the while loop:
                # Read lines from thisSample and change in thisOut. Hence for a new iteration I need to sync the two.
                thisSample = thisOut

                rleRaw = rle(thisSample$smooth)
                # This takes the indeces of the FALSE chains and adds 1. This should give you the next segment which is TRUE.
                # Two challenges:
                # 1) Last segment always FALSE (see above), hence removal of the last number as this would indicate to a segment outside the df.
                # 2) If it starts with a TRUE segment, this would not be found when looking at the FALSE chains. Hence, adding index 1 manually if chain starts with TRUE.
                indRaw = cumsum(rleRaw$lengths)[ ! rleRaw$values ] + 1
                indRaw = indRaw[ -length(indRaw) ]
                if( rleRaw$values[1] ) { indRaw = c(1, indRaw) }

                # loop over start indices of TRUE chains.
                for(i in indRaw) {
                    # detect length of segments to smooth. add 1 as the last segment has a FALSE value in it but still belongs to this chain.
                    endOfStreak = i + rle(thisSample$smooth[i:nrow(thisSample)])$lengths[1]
                    # extract reads
                    dfMerge = thisSample[i:endOfStreak,]

                    # too stupid to make this work with data.table
                    newElement = as.data.frame( dfMerge[1,] )
                    # Get new end and check first wether valid number.
                    newEnd = dfMerge[nrow(dfMerge),][[colNameEnd]]
                    if(! is.null(newEnd)) {
                        newElement[[colNameEnd]] = newEnd
                    } else {
                        stop("New end coordinate is null. Supplied correct column name?")
                    }
                    ## Column "segVal" will be dealt with in a minute. Column "diffs" later when running again idSmoothingTargets.

                    # Merge cn specifically by taking the length of the elements into consideration
                    widthWeights = dfMerge[[colNameEnd]] - dfMerge[[colNameStart]]
                    newElement[[colNameMerge[1]]] = stats::weighted.mean(dfMerge[[colNameMerge[1]]], widthWeights)
                    # Replace all to merge segments with the new merged segment. Later delete duplicated.
                    thisOut[i:endOfStreak,] = newElement
                }

                # as we have replaced all segments with the new mean segment, we need to remove the duplicates
                thisOut = thisOut[ ! duplicated(thisOut), ]
                # again detect segments which needs smoothing
                thisOut = idSmoothingTargets(thisOut, WIGGLE, colNameSegVal = colNameMerge[[1]], colNameChr = colNameChr,
                                             IGNOREDELS = IGNOREDELS)
                stillSmoothing = sum(thisOut$smooth)
            }

            # after smoothing is finished, change name of cohort
            thisOut$smooth = NULL
            thisOut$diffs = NULL
            return( thisOut )
        }
        doParallel::stopImplicitCluster()
    } else {
        #stop("no single thread method for smoothing")
        lSmooth <- lapply(lRaw,FUN = function(x){

            thisOut = x
            stillSmoothing = sum(thisOut$smooth)
            while( stillSmoothing > 0 ) {
                # For the while loop:
                # Read lines from thisSample and change in thisOut. Hence for a new iteration I need to sync the two.
                thisSample = thisOut

                rleRaw = rle(thisSample$smooth)
                # This takes the indeces of the FALSE chains and adds 1. This should give you the next segment which is TRUE.
                # Two challenges:
                # 1) Last segment always FALSE (see above), hence removal of the last number as this would indicate to a segment outside the df.
                # 2) If it starts with a TRUE segment, this would not be found when looking at the FALSE chains. Hence, adding index 1 manually if chain starts with TRUE.
                indRaw = cumsum(rleRaw$lengths)[ ! rleRaw$values ] + 1
                indRaw = indRaw[ -length(indRaw) ]
                if( rleRaw$values[1] ) { indRaw = c(1, indRaw) }

                # loop over start indices of TRUE chains.
                for(i in indRaw) {
                    # detect length of segments to smooth. add 1 as the last segment has a FALSE value in it but still belongs to this chain.
                    endOfStreak = i + rle(thisSample$smooth[i:nrow(thisSample)])$lengths[1]
                    # extract reads
                    dfMerge = thisSample[i:endOfStreak,]

                    # too stupid to make this work with data.table
                    newElement = as.data.frame( dfMerge[1,] )
                    # Get new end and check first wether valid number.
                    newEnd = dfMerge[nrow(dfMerge),][[colNameEnd]]
                    if(! is.null(newEnd)) {
                        newElement[[colNameEnd]] = newEnd
                    } else {
                        stop("New end coordinate is null. Supplied correct column name?")
                    }
                    ## Column "segVal" will be dealt with in a minute. Column "diffs" later when running again idSmoothingTargets.

                    # Merge cn specifically by taking the length of the elements into consideration
                    widthWeights = dfMerge[[colNameEnd]] - dfMerge[[colNameStart]]
                    newElement[[colNameMerge[1]]] = stats::weighted.mean(dfMerge[[colNameMerge[1]]], widthWeights)
                    # Replace all to merge segments with the new merged segment. Later delete duplicated.
                    thisOut[i:endOfStreak,] = newElement
                }

                # as we have replaced all segments with the new mean segment, we need to remove the duplicates
                thisOut = thisOut[ ! duplicated(thisOut), ]
                # again detect segments which needs smoothing
                thisOut = idSmoothingTargets(thisOut, WIGGLE, colNameSegVal = colNameMerge[[1]], colNameChr = colNameChr,
                                             IGNOREDELS = IGNOREDELS)
                stillSmoothing = sum(thisOut$smooth)
            }

            # after smoothing is finished, change name of cohort
            thisOut$smooth = NULL
            thisOut$diffs = NULL
            return( thisOut )
        })
    }
    if( isTRUE(asDf) ) {
        dfSmooth = setDT( rbindlist( lSmooth ) )
        return( dfSmooth )
    } else {
        return( lSmooth )
    }

}
