# clinPredictionPlatinum

data(SIGobj.drews)
data(SIGobj.mac)
data(CNobj.drews)

out <- rep("Predicted resistant",times=length(rownames(SIGobj.drews@samplefeatData)))
names(out) <- rownames(SIGobj.drews@samplefeatData)

test_that("test no sig", {
    expect_error(clinPredictionPlatinum(object = CNobj.drews))
})

test_that("test method", {
    expect_error(clinPredictionPlatinum(object = SIGobj.mac))
})

test_that("test correct outcome", {
    expect_equal(clinPredictionPlatinum(object = SIGobj.drews),out)
})
