# clinPredictionDenovo

data(SIGobj.drews)
data(SIGobj.mac)
data(CNobj.drews)

test_that("function depreciated", {
  expect_error(clinPredictionDenovo())
})

test_that("function depreciated", {
    expect_error(clinPredictionDenovo(object = CNobj.drews))
})

test_that("function depreciated", {
    expect_error(clinPredictionDenovo(object = SIGobj.mac))
})

test_that("function depreciated", {
    expect_error(clinPredictionDenovo(object = SIGobj.drews))
})
