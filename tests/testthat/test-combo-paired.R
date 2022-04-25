# This tests countPairedComboBarcodes().
# library(testthat); library(screenCounter); source("setup.R"); source("test-combo-paired.R")

set.seed(100000)
library(Biostrings)
nbarcodes1 <- 20
POOL1 <- vapply(rep(10, nbarcodes1), GENERATE_RANDOM_SEQ, FUN.VALUE="")
nbarcodes2 <- 50
POOL2 <- vapply(rep(15, nbarcodes2), GENERATE_RANDOM_SEQ, FUN.VALUE="")

DUMPSEQ <- function(fname, n, barcode.fmt, pool, strandFUN = identity) {
    i <- sample(length(pool), n, replace=TRUE)
    barcodes <- sprintf(barcode.fmt, pool[i])
    names(barcodes) <- seq_len(n)
    ADD_FLANKS(barcodes, fname, strandFUN=strandFUN)
    i
}

CHECK_OUTPUT <- function(out, ref1, ref2, count, N) {
    expect_identical(out$combinations[,1], ref1)
    expect_identical(out$combinations[,2], ref2)
    expect_identical(out$counts, count)
    expect_identical(metadata(out)$npairs, as.integer(N))
}

barcode.fmt1 <- "ACGT%sACGT"
template.fmt1 <- sub("%s", strrep("-", 10), barcode.fmt1)
barcode.fmt2 <- "AAAA%sTTTT"
template.fmt2 <- sub("%s", strrep("-", 15), barcode.fmt2)

test_that("countPairedComboBarcodes works as expected in basic mode", {
    tmp1 <- tempfile()
    tmp2 <- tempfile()

    for (N in c(1, 10, 100, 1000)) {
        i1 <- DUMPSEQ(tmp1, N, barcode.fmt1, POOL1)
        i2 <- DUMPSEQ(tmp2, N, barcode.fmt2, POOL2)
        out <- countPairedComboBarcodes(c(tmp1, tmp2), template=c(template.fmt1, template.fmt2), choices=list(POOL1, POOL2), indices=TRUE)

        ref <- REF_COMBINE(i1, i2, POOL1, POOL2)
        CHECK_OUTPUT(out, ref[[1]], ref[[2]], ref[[3]], N)
    }
})

test_that("countPairedComboBarcodes works as expected with strands", {
    N <- 1000
    tmp1 <- tempfile()
    tmp2 <- tempfile()

    # Original and reverse.
    i1 <- DUMPSEQ(tmp1, N, barcode.fmt1, POOL1)
    i2 <- DUMPSEQ(tmp2, N, barcode.fmt2, POOL2, CHOOSE_STRAND_FUN("reverse"))
    out <- countPairedComboBarcodes(c(tmp1, tmp2), template=c(template.fmt1, template.fmt2), choices=list(POOL1, POOL2), strand=c("original", "reverse"), indices=TRUE)

    ref <- REF_COMBINE(i1, i2, POOL1, POOL2)
    CHECK_OUTPUT(out, ref[[1]], ref[[2]], ref[[3]], N)

    # Reverse and original
    i1 <- DUMPSEQ(tmp1, N, barcode.fmt1, POOL1, CHOOSE_STRAND_FUN("reverse"))
    i2 <- DUMPSEQ(tmp2, N, barcode.fmt2, POOL2)
    out <- countPairedComboBarcodes(c(tmp1, tmp2), template=c(template.fmt1, template.fmt2), choices=list(POOL1, POOL2), strand=c("reverse",  "original"), indices=TRUE)

    ref <- REF_COMBINE(i1, i2, POOL1, POOL2)
    CHECK_OUTPUT(out, ref[[1]], ref[[2]], ref[[3]], N)
})

test_that("countPairedComboBarcodes works with reporting sequences directly", {
    N <- 100
    tmp1 <- tempfile()
    tmp2 <- tempfile()

    i1 <- DUMPSEQ(tmp1, N, barcode.fmt1, POOL1)
    i2 <- DUMPSEQ(tmp2, N, barcode.fmt2, POOL2)

    out1 <- countPairedComboBarcodes(c(tmp1, tmp2), template=c(template.fmt1, template.fmt2), choices=list(POOL1, POOL2))
    out2 <- countPairedComboBarcodes(c(tmp1, tmp2), template=c(template.fmt1, template.fmt2), choices=list(POOL1, POOL2), indices=TRUE)

    expect_identical(out1$combinations[,1], POOL1[out2$combinations[,1]])
    expect_identical(out1$combinations[,2], POOL2[out2$combinations[,2]])
})

test_that("matrixOfPairedComboBarcodes handles names correctly", {
    N <- 100
    tmp1a <- tempfile()
    tmp1b <- tempfile()
    i1 <- DUMPSEQ(tmp1a, N, barcode.fmt1, POOL1)
    i2 <- DUMPSEQ(tmp1b, N, barcode.fmt2, POOL2)

    # Generating a second file.
    N <- 200
    tmp2a <- tempfile()
    tmp2b <- tempfile()
    i1 <- DUMPSEQ(tmp2a, N, barcode.fmt1, POOL1)
    i2 <- DUMPSEQ(tmp2b, N, barcode.fmt2, POOL2)

    # File names match up.
    files <- list(c(tmp1a, tmp1b), c(tmp2a, tmp2b))
    se <- matrixOfPairedComboBarcodes(files, template=c(template.fmt1, template.fmt2), list(first=POOL1, second=POOL2))
    expect_identical(colnames(se), basename(c(tmp1a, tmp2a)))

    # Works when pool names are not specified.
    se <- matrixOfPairedComboBarcodes(files, template=c(template.fmt1, template.fmt2), list(POOL1, POOL2))
    expect_identical(colnames(rowData(se)), c("first", "second"))
}) 
