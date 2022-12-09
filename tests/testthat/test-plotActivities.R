## plotActivities
# signature activities

test_that("null input", {
  expect_error(plotActivities())
})

test_that("class input error", {
    expect_error(plotActivities("data"))
})
