# checkbinned

non.binned <- data.frame(chromosome=c(1,1),start=c(1,1001),end=c(1000,2500))
binned <- data.frame(chromosome=c(1,1),start=c(1,1001),end=c(1000,2000))

test_that("test no data", {
    expect_error(checkbinned())
})

test_that("test binned", {
  expect_equal(checkbinned(binned),TRUE)
})

test_that("test non-binned", {
    expect_equal(checkbinned(non.binned),FALSE)
})
