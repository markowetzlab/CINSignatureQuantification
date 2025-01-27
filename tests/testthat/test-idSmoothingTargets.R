## idSmoothingTargets

# idSmoothingTargets = function(dfAllSegs, WIGGLE, colNameSegVal, colNameChr, IGNOREDELS = TRUE) {
#
#     ### Check column name
#     testSegVal = dfAllSegs[[colNameSegVal]][1]
#     testChr = dfAllSegs[[colNameChr]][1]
#
#     # Quick sanity checks
#     if(! is.numeric(testSegVal)) { stop("Segment Value column has no numeric value in it. Supplied correct column name? Forgot conversion?")}
#     if(is.null(testSegVal)) { stop("Chromosome column has no numeric value in it. Supplied correct column name?")}
#
#     # Take differences to segment down below
#     dfAllSegs$diffs = c( abs( dfAllSegs[[colNameSegVal]][1:(nrow(dfAllSegs)-1)] - dfAllSegs[[colNameSegVal]][2:nrow(dfAllSegs)] ), WIGGLE+1)
#     # Set TRUE if difference to next segment is smaller than the user supplied cutoff
#     dfAllSegs$smooth = dfAllSegs$diffs <= WIGGLE
#     # Set all segments which are last in a chromosome to FALSE. This also prevents leaking to other samples and cohorts.
#     dfAllSegs$smooth[ cumsum( rle(as.character(dfAllSegs[[colNameChr]]))$lengths ) ] = FALSE
#
#     # Ignore deletions if wished
#     if(IGNOREDELS) { dfAllSegs$smooth[ dfAllSegs[[colNameSegVal]] == 0 ] = FALSE }
#
#     return( dfAllSegs )
# }

data("TCGA_478_Samples_SNP6_GOLD")
dfAllSegs = as.data.frame(TCGA_478_Samples_SNP6_GOLD)
WIGGLE = 0.1
colNameSegVal = "segVal"
colNameChr = "chromosome"
IGNOREDELS = FALSE

dfAllSegsTest <- dfAllSegs
dfAllSegsTest$diffs = c( abs( dfAllSegsTest[[colNameSegVal]][1:(nrow(dfAllSegsTest)-1)] - dfAllSegsTest[[colNameSegVal]][2:nrow(dfAllSegsTest)] ), WIGGLE+1)
dfAllSegsTest$smooth = dfAllSegsTest$diffs <= WIGGLE
dfAllSegsTest$smooth[ cumsum( rle(as.character(dfAllSegsTest[[colNameChr]]))$lengths ) ] = FALSE

test_that("empty input", {
  expect_error(idSmoothingTargets())
})

test_that("default outcome", {
    expect_equal(idSmoothingTargets(dfAllSegs,
                                    WIGGLE = 0.1,
                                    colNameSegVal = "segVal",
                                    colNameChr = "chromosome",
                                    IGNOREDELS = IGNOREDELS),
    dfAllSegsTest)
})

test_that("wrong segVal type", {
    expect_error(idSmoothingTargets(dfAllSegs,
                                    WIGGLE = 0.1,
                                    colNameSegVal = "chromosome",
                                    colNameChr = "chromosome",
                                    IGNOREDELS = IGNOREDELS))
})

test_that("null segVal type", {
    expect_error(idSmoothingTargets(dfAllSegs,
                                    WIGGLE = 0.1,
                                    colNameSegVal = "s",
                                    colNameChr = "chromosome",
                                    IGNOREDELS = IGNOREDELS))
})

dfAllSegsTest$smooth[ dfAllSegsTest[[colNameSegVal]] == 0 ] = FALSE

test_that("IGNORE deletions", {
    expect_equal(idSmoothingTargets(dfAllSegs,
                                    WIGGLE = 0.1,
                                    colNameSegVal = "segVal",
                                    colNameChr = "chromosome",
                                    IGNOREDELS = TRUE),
    dfAllSegsTest)
})
