test_that("multiplication works", {
  expect_equal(checkSegValRounding(c(1.8,4.0,5.3)),FALSE)
  expect_equal(checkSegValRounding(c(1,0,5)),TRUE)
})
