# scaleByModel

newModel <- matrix(runif(30),ncol = 3,
                   dimnames = list(paste0("sample",seq.int(1,10,1)),c("S1","S2","S3")))
newModel <- scale(newModel)
newModel <- list(mean = attr(newModel, "scaled:center"),
                scale = attr(newModel, "scaled:scale"))

newmNorm <- matrix(runif(30),ncol = 3,
          dimnames = list(paste0("sample",seq.int(1,10,1)),c("S1","S2","S3")))

testOut <- sweep(newmNorm, 2, newModel$mean, FUN = '-')
testOut <- sweep(testOut,2, newModel$scale, FUN = "/")

test_that("test no data", {
    expect_error(scaleByModel())
})

test_that("test one input", {
    expect_error(scaleByModel(newmNorm))
})

test_that("test one input", {
    expect_error(scaleByModel(newModel))
})

test_that("test matrix return", {
    expect_true(is.matrix(scaleByModel(newmNorm,newModel)))
})

test_that("test return", {
    expect_equal(scaleByModel(newmNorm,newModel),testOut)
})
