# getSamples

data(CNobj.drews)
t <- rownames(CNobj.drews@samplefeatData)

test_that("return type", {
  expect_type(getSamples(CNobj.drews),"character")
})

test_that("test identical", {
    expect_equal(getSamples(CNobj.drews),t)
})
