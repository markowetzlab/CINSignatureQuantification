# clinPredictionDenovo

data(SIGobj.drews)
data(SIGobj.mac)
data(CNobj.drews)

test_that("empty input", {
    expect_error(clinPredictionDenovo())
})

test_that("CNQuat error", {
    expect_error(clinPredictionDenovo(object = CNobj.drews))
})

test_that("function depreciated mac", {
    expect_message(clinPredictionDenovo(object = SIGobj.mac))
})

test_that("function depreciated drews", {
    expect_message(clinPredictionDenovo(object = SIGobj.drews))
})
