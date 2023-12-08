# getActivities

data(SIGobj.drews)
data(SIGobj.mac)

td <- SIGobj.drews@activities$thresholdAct2
tm <- SIGobj.mac@activities$thresholdAct2

test_that("test incorrect type drews", {
    expect_error(getActivities(SIGobj.drews,type = "unsupported"))
})

test_that("test unsupported type mac", {
    expect_error(getActivities(SIGobj.mac,type = "norm"))
    expect_error(getActivities(SIGobj.mac,type = "scaled"))
})

test_that("test return equal default", {
    expect_equal(getActivities(SIGobj.drews),td)
    expect_equal(getActivities(SIGobj.mac),tm)
})

td <- SIGobj.drews@activities$rawAct0
tm <- SIGobj.mac@activities$rawAct0

test_that("test return equal raw", {
    expect_equal(getActivities(SIGobj.drews,type = "raw"),td)
    expect_equal(getActivities(SIGobj.mac,type = "raw"),tm)
})

td <- SIGobj.drews@activities$normAct1

test_that("test return equal norm", {
    expect_equal(getActivities(SIGobj.drews,type = "norm"),td)
})

td <- SIGobj.drews@activities$scaledAct3

test_that("test return equal scaled", {
    expect_equal(getActivities(SIGobj.drews,type = "scaled"),td)
})
