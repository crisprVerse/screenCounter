# This tests countComboBarcodes().
# library(testthat); library(screenCounter); source("setup.R"); source("test-combo-single.R")

library(Biostrings)
POOL1 <- vapply(rep(10, 20), GENERATE_RANDOM_SEQ, FUN.VALUE="")
POOL2 <- vapply(rep(8, 10), GENERATE_RANDOM_SEQ, FUN.VALUE="") 

barcode.fmt <- "ACGT%sACGT%sACGT"
template <- sprintf(barcode.fmt, strrep("N", nchar(POOL1[1])), strrep("N", nchar(POOL2[1])))

STICKER <- function(barcodes, fname, out, ..., strandFUN=identity) {
    ADD_FLANKS(barcodes, fname, strandFUN=strandFUN)
    out2 <- countComboBarcodes(fname, template, ..., indices=TRUE)
    expect_identical(out, out2)
}

CHECK_OUTPUT <- function(out, ref1, ref2, count, N) {
    expect_identical(out$combinations[,1], ref1)
    expect_identical(out$combinations[,2], ref2)
    expect_identical(out$counts, count)
    expect_identical(metadata(out)$nreads, as.integer(N))
}

test_that("countComboBarcodes works as expected in basic mode", {
    for (N in c(1, 10, 100, 1000)) {
        i1 <- sample(length(POOL1), N, replace=TRUE)
        i2 <- sample(length(POOL2), N, replace=TRUE)
        barcodes <- sprintf(barcode.fmt, POOL1[i1], POOL2[i2])
        names(barcodes) <- seq_len(N)

        tmp <- tempfile(fileext=".fastq")
        writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

        out <- countComboBarcodes(tmp, template, list(POOL1, POOL2), indices=TRUE)
        ref <- REF_COMBINE(i1, i2, POOL1, POOL2)
        CHECK_OUTPUT(out, ref[[1]], ref[[2]], ref[[3]], N)

        # Same results when you stick a bunch of random crap to the start and end.
        STICKER(barcodes, tmp, out, choices=list(POOL1, POOL2))

        # Should be the same, just checking that it works.
        best <- countComboBarcodes(tmp, template, list(POOL1, POOL2), find.best=TRUE, indices=TRUE)
        expect_identical(out, best)
    }
})

test_that("countComboBarcodes works as expected with strands", {
    N <- 1000
    for (strand in c("original", "reverse", "both")) {
        i1 <- sample(length(POOL1), N, replace=TRUE)
        i2 <- sample(length(POOL2), N, replace=TRUE)
        barcodes <- sprintf(barcode.fmt, POOL1[i1], POOL2[i2])
        names(barcodes) <- seq_len(N)

        B <- DNAStringSet(barcodes)
        strandFUN <- CHOOSE_STRAND_FUN(strand)
        tmp <- tempfile(fileext=".fastq")
        writeXStringSet(strandFUN(B), filepath=tmp, format="fastq")

        out <- countComboBarcodes(tmp, template, list(POOL1, POOL2), strand=strand, indices=TRUE)
        ref <- REF_COMBINE(i1, i2, POOL1, POOL2)
        CHECK_OUTPUT(out, ref[[1]], ref[[2]], ref[[3]], N)

        # Same results when you stick a bunch of random crap to the start and end.
        STICKER(barcodes, tmp, out, choices=list(POOL1, POOL2), strand=strand, strandFUN=strandFUN)
    }
})

test_that("countComboBarcodes works with reporting sequences directly", {
    N <- 100
    i1 <- sample(length(POOL1), N, replace=TRUE)
    i2 <- sample(length(POOL2), N, replace=TRUE)
    barcodes <- sprintf(barcode.fmt, POOL1[i1], POOL2[i2])
    names(barcodes) <- seq_len(N)

    B <- DNAStringSet(barcodes)
    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(B, filepath=tmp, format="fastq")

    out1 <- countComboBarcodes(tmp, template, list(POOL1, POOL2), indices=FALSE)
    out2 <- countComboBarcodes(tmp, template, list(POOL1, POOL2), indices=TRUE)
    expect_identical(out1$combinations[,1], POOL1[out2$combinations[,1]])
    expect_identical(out1$combinations[,2], POOL2[out2$combinations[,2]])
})

test_that("matrixOfComboBarcodes handles names correctly", {
    N <- 100
    i1 <- sample(length(POOL1), N, replace=TRUE)
    i2 <- sample(length(POOL2), N, replace=TRUE)
    barcodes <- sprintf(barcode.fmt, POOL1[i1], POOL2[i2])
    names(barcodes) <- seq_len(N)

    B <- DNAStringSet(barcodes)
    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(B, filepath=tmp, format="fastq")

    # Generating a second file.
    N <- 200
    i1 <- sample(length(POOL1), N, replace=TRUE)
    i2 <- sample(length(POOL2), N, replace=TRUE)
    barcodes <- sprintf(barcode.fmt, POOL1[i1], POOL2[i2])
    names(barcodes) <- seq_len(N)

    B <- DNAStringSet(barcodes)
    tmp2 <- tempfile(fileext=".fastq")
    writeXStringSet(B, filepath=tmp2, format="fastq")

    # File names match up.
    se <- matrixOfComboBarcodes(c(tmp, tmp2), template, list(first=POOL1, second=POOL2))
    expect_identical(colnames(se), basename(c(tmp, tmp2)))

    # Works when pool names are not specified.
    se <- matrixOfComboBarcodes(tmp, template, list(POOL1, POOL2))
    expect_identical(colnames(rowData(se)), c("first", "second"))
}) 
