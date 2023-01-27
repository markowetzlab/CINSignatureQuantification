# subsetSigActivities

# . <- function(x,s){
#     subSigActivities <- lapply(x, FUN = function(y){
#         y <- y[which(rownames(y) %in% s),]
#     })
#     subSigActivities
# }

data(SIGobj.drews)
x <- SIGobj.drews@activities
s <- rownames(x)[1:10]

xs <- lapply(x, FUN = function(y){
            y <- y[which(rownames(y) %in% s),]
        })

test_that("empty input", {
    expect_error(subsetSigActivities())
})

test_that("correct output", {
    expect_equal(subsetSigActivities(x,s),xs)
})
