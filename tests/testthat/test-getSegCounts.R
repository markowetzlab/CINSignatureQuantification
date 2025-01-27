# getSegCounts

data(CNobj.drews)
cn <- CNobj.drews@segments
counts <- lapply(cn,nrow)
counts <- unlist(counts)

test_that("test no input", {
  expect_error(getSegCounts())
})

test_that("test output", {
    expect_equal(getSegCounts(cn),counts)
})

test_that("test return type", {
    expect_type(getSegCounts(cn),"integer")
})

test_that("test named", {
    expect_named(getSegCounts(cn))
})
