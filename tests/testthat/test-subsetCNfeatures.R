# subsetCNfeatures

data(SIGobj.drews)
x <- SIGobj.drews@featData
s <- unique(x$segsize$ID)[1:10]

xs <- lapply(x, FUN = function(y){
    y <- y[y$ID %in% s,]
})

test_that("empty input", {
    expect_error(subsetCNfeatures())
})

test_that("correct output SIGQuant.drews", {
    expect_equal(subsetCNfeatures(x,s),xs)
})

data(CNobj.drews)
x <- CNobj.drews@featData
s <- unique(x$segsize$ID)[1:10]

xs <- lapply(x, FUN = function(y){
    y <- y[y$ID %in% s,]
})

test_that("correct output CNQuant.drews", {
    expect_equal(subsetCNfeatures(x,s),xs)
})

data(CNobj.mac)
x <- CNobj.mac@featData
s <- unique(x$segsize$ID)[1:10]

xs <- lapply(x, FUN = function(y){
    y <- y[y$ID %in% s,]
})

test_that("correct output CNQuant.mac", {
    expect_equal(subsetCNfeatures(x,s),xs)
})

data(SIGobj.mac)
x <- SIGobj.mac@featData
s <- unique(x$segsize$ID)[1:10]

xs <- lapply(x, FUN = function(y){
    y <- y[y$ID %in% s,]
})

test_that("correct output SIGobj.mac", {
    expect_equal(subsetCNfeatures(x,s),xs)
})
