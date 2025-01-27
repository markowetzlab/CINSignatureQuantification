## CheckChromosomeFormat
# Set up chromosomes
chrs <- c(seq.int(1:22),c("X"))
chrs_chr <- c(paste0("chr",c(seq.int(1:22),c("X"))))
chrs_num <- c(seq.int(1:23))

test_that("check no chromosome formatting", {
    expect_equal(checkChromosomeFormat(chrs),chrs)
})

test_that("check chromosome prefix", {
    expect_message(expect_equal(checkChromosomeFormat(chrs_chr),chrs))
})

test_that("check numerical non-autosomes", {
    expect_message(expect_equal(checkChromosomeFormat(chrs_num),chrs))
})

test_that("check null input", {
    expect_error(checkChromosomeFormat())
})
