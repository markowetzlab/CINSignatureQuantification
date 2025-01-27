# subsetfeatFitting

data(SIGobj.drews)
x <- SIGobj.drews@featFitting
s <- rownames(x$sampleByComponent)[1:10]

xs <- x
sxc <- x$sampleByComponent
sxc <- sxc[which(rownames(sxc) %in% s),]
xs$sampleByComponent <- sxc

test_that("empty input", {
    expect_error(subsetfeatFitting())
})

test_that("correct output SIGQuant.drews", {
    expect_equal(subsetfeatFitting(x,s),xs)
})
