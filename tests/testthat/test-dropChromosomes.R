## dropChromosomes
# Set up chromosomes
chrs <- c(seq.int(1:22),c("X"))
d <- data.frame(chromosome=c("1","2","Y","mt"),
                start=c(1,1,1,1),
                end=c(2,2,2,2),
                segVal=c(1,2,3,4),
                sample=c("s1","s1","s1","s1"))
t <- d[d$chromosome %in% c("1","2"),]
t$chromosome <- factor(t$chromosome,levels = chrs)

test_that("check drop Chromosomes", {
  expect_message(expect_equal(dropChromosomes(d),t))
})

test_that("check null input", {
  expect_error(dropChromosomes())
})
