# getSampleByComponent

data(CNobj.drews)
t <- CNobj.drews@featFitting$sampleByComponent

test_that("test identical", {
    expect_equal(getSampleByComponent(CNobj.drews),t)
})
