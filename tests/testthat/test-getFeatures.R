# getFeatures

data(CNobj.drews)
CNobj.drews.null <- CNobj.drews
CNobj.drews.null@featData <- list()

test_that("test no data", {
    expect_error(getFeatures(CNobj.drews.null))
})

f <- CNobj.drews@featData
fs <- CNobj.drews@featData["segsize"]

test_that("test return all", {
    expect_equal(getFeatures(CNobj.drews),f)
})

test_that("test return subset", {
    expect_equal(getFeatures(CNobj.drews,feat = "segsize"),fs)
})

test_that("test return subset", {
    expect_error(getFeatures(CNobj.drews,feat = "unspecified"))
})
