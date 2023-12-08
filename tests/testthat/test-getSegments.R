# getSegments

data(CNobj.drews)
segTable <- do.call(rbind,CNobj.drews@segments)
t <- data.frame(segTable,row.names = NULL)

test_that("test identical", {
    expect_equal(getSegments(CNobj.drews),t)
})

test_that("test return type", {
    expect_type(getSegments(CNobj.drews),"list")
})
