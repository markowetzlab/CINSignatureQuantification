# removeQuietSamples

dtSmooth <- do.call(rbind,lapply(SIGobj.drews@segments,
                                 FUN = function(x) return(x)))

test_that("test empty", {
  expect_error(removeQuietSamples())
})

test_that("test default", {
    expect_equal(removeQuietSamples(dtSmooth),dtSmooth)
})

test_that("test DCIN input string", {
    expect_error(removeQuietSamples(dtSmooth,DCIN = "1"))
})

test_that("test DCIN input float", {
    expect_error(removeQuietSamples(dtSmooth,DCIN = 2.4))
})

test_that("test CDIN input -Ve", {
    expect_error(removeQuietSamples(dtSmooth,DCIN = -2))
})

quietSamples = c("TCGA-BT-A20P","TCGA-BT-A20Q","TCGA-BT-A3PH","TCGA-CF-A27C","TCGA-CF-A3MF")
dtSmooth.1 = dtSmooth[ ! dtSmooth$sample %in% quietSamples, ]
dtSmooth.1$sample = factor(dtSmooth.1$sample)

test_that("test DCIN", {
    expect_equal(removeQuietSamples(dtSmooth,DCIN = 100),dtSmooth.1)
})

