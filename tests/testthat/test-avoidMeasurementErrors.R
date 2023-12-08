# avoidMeasurementErrors

data("CNobj.drews")
cnTable <- do.call(rbind,CNobj.drews@segments)
cnTableT <- cnTable

# add zero errors
cnTable$segVal[which(cnTable$segVal < 0)] <- round(runif(8,min = -10,max = -0.1),
                                                   digits = 2)
# Non-factor samples
cnTable$sample <- as.character(cnTable$sample)

test_that("test no input", {
  expect_error(avoidMeasurementErrors())
})

test_that("test correct output", {
    expect_equal(avoidMeasurementErrors(cnTable),cnTableT)
})

