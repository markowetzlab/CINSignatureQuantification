# getExperiment

data(CNobj.drews)
t <- CNobj.drews@ExpData

test_that("test identical", {
    expect_equal(getExperiment(CNobj.drews),t)
})
