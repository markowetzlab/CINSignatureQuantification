# getSamplefeatures

data(CNobj.drews)
t <- as.data.frame(CNobj.drews@samplefeatData)

test_that("return type", {
    expect_type(getSamplefeatures(CNobj.drews),"list")
})

test_that("test identical", {
    expect_equal(getSamplefeatures(CNobj.drews),t)
})
