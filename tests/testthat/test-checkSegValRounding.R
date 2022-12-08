## checkSegValRounding
test_that("check rounding false", {
    expect_equal(checkSegValRounding(c(1.8,4.0,5.3)),FALSE)
})

test_that("check rounding true", {
    expect_equal(checkSegValRounding(c(1,0,5)),TRUE)
})
