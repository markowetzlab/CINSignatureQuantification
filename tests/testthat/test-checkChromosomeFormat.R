test_that("checkChromosomeFormat", {
    chrs <- c(seq.int(1:22),c("X"))
    chrs_chr <- c(paste0("chr",c(seq.int(1:22),c("X"))))
    chrs_num <- c(seq.int(1:23))

    expect_equal(checkChromosomeFormat(chrs),chrs)
    expect_equal(checkChromosomeFormat(chrs_chr),chrs)
    expect_equal(checkChromosomeFormat(chrs_num),chrs)
})
