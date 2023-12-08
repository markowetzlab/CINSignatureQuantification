# getPloidyfeat

data("CNobj.drews")
x <- CNobj.drews@segments

featploidy <- unlist(lapply(x,FUN = function(y){
    segLen<-(as.numeric(y$end)-as.numeric(y$start))
    ploidy<-sum((segLen/sum(segLen))*as.numeric(y$segVal))
}))
featploidy <- round(featploidy,digits = 3)

test_that("test no input", {
    expect_error(getPloidyfeat())
})

test_that("test output", {
    expect_equal(getPloidyfeat(x),featploidy)
})

test_that("test return type", {
    expect_type(getPloidyfeat(x),"double")
})

test_that("test named", {
    expect_named(getPloidyfeat(x))
})
