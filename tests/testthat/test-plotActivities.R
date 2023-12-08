## plotActivities
# signature activities

data("SIGobj.drews")
data("SIGobj.mac")

test_that("null input", {
  expect_error(plotActivities())
})

test_that("class input error", {
    expect_error(plotActivities("data"))
})

test_that("unknown method", {
    expect_error(plotActivities(object = SIGobj.drews,type = "unknown"))
})

test_that("drews method threshold", {
    expect_type(plotActivities(object = SIGobj.drews,type = "threshold"),"list")
    expect_no_condition(plotActivities(object = SIGobj.drews,type = "threshold"))
})

test_that("mac method threshold", {
    expect_type(plotActivities(object = SIGobj.mac,type = "threshold"),"list")
    expect_no_condition(plotActivities(object = SIGobj.mac,type = "threshold"))
})

test_that("drews method raw", {
    expect_type(plotActivities(object = SIGobj.drews,type = "raw"),"list")
    expect_no_condition(plotActivities(object = SIGobj.drews,type = "raw"))
})

test_that("mac method raw", {
    expect_type(plotActivities(object = SIGobj.mac,type = "raw"),"list")
    expect_no_condition(plotActivities(object = SIGobj.mac,type = "raw"))
})

test_that("drews method norm", {
    expect_type(plotActivities(object = SIGobj.drews,type = "norm"),"list")
    expect_no_condition(plotActivities(object = SIGobj.drews,type = "norm"))
})

test_that("mac method norm", {
    expect_error(plotActivities(object = SIGobj.mac,type = "norm"))
})

test_that("drews method scaled", {
    expect_type(plotActivities(object = SIGobj.drews,type = "scaled"),"list")
    expect_no_condition(plotActivities(object = SIGobj.drews,type = "scaled"))
})

test_that("mac method scaled", {
    expect_error(plotActivities(object = SIGobj.mac,type = "scaled"))
})


test_that("custom colours", {
    expect_error(plotActivities(object = SIGobj.drews,cols = c("1","2")))
})
