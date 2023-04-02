# This tests countRandomBarcodes().
# library(testthat); library(screenCounter); source("setup.R"); source("test-random.R")

library(Biostrings)
len <- 6
barcode.fmt <- "AAAAACGT%sACGTGGGG" # longer template to avoid matches to random flanks; also non-palindromic, to avoid ambiguity when best = TRUE.
template <- sprintf(barcode.fmt, strrep("N", len))

STICKER <- function(barcodes, fname, out, ..., strandFUN=identity) {
    ADD_FLANKS(barcodes, fname, strandFUN=strandFUN)
    out2 <- countRandomBarcodes(fname, template=template, ...)
    expect_identical(out, out2)

    # Nothing at the start, flanking on the right.
    ADD_FLANKS(barcodes, fname, nleft=0, strandFUN=strandFUN)
    out2 <- countRandomBarcodes(fname, template=template, ...)
    expect_identical(out, out2)

    # Nothing at the end, flanking on the left.
    ADD_FLANKS(barcodes, fname, nright=0, strandFUN=strandFUN)
    out2 <- countRandomBarcodes(fname, template=template, ...)
    expect_identical(out, out2)
}

test_that("countRandomBarcodes works as expected in basic mode", {
    for (N in c(1, 10, 100, 1000)) {
        randomized <- vapply(rep(len, N), GENERATE_RANDOM_SEQ, FUN.VALUE="")
        barcodes <- sprintf(barcode.fmt, randomized)
        names(barcodes) <- seq_len(N)

        tmp <- tempfile(fileext=".fastq")
        writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

        out <- countRandomBarcodes(tmp, template=template)
        expect_identical(as.integer(N), metadata(out)$nreads) 

        tab <- table(randomized)
        o <- order(names(tab))
        expect_identical(names(tab)[o], out$sequences)
        expect_identical(as.integer(unname(tab))[o], out$counts)

        # Same results when you stick a bunch of random crap to the start and end.
        STICKER(barcodes, tmp, out)

        # Same results with find.best=TRUE.
        best <- countRandomBarcodes(tmp, template=template, find.best=TRUE)
        expect_identical(out, best)
    }
})

test_that("countRandomBarcodes works as expected with substitutions", {
    barcodes <- c(
        "AAAAACGTGGGGGGACGTGGGG",
        "AAATACGTGGGGGCACGTGGGG",
        "AAATACGTGGGGCCACGTGGGC" # not matched, two substitutions.
    )
    names(barcodes) <- seq_along(barcodes)

    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

    out <- countRandomBarcodes(tmp, template=template, substitutions=1)
    expect_identical(out$sequences, c("GGGGGC", "GGGGGG"))
    expect_identical(out$counts, c(1L, 1L))

    out <- countRandomBarcodes(tmp, template=template, substitutions=2)
    expect_identical(out$sequences, c("GGGGCC", "GGGGGC", "GGGGGG"))
    expect_identical(out$counts, c(1L, 1L, 1L))

    STICKER(barcodes, tmp, out, substitutions=2)
})

test_that("countRandomBarcodes works as expected with strands", {
    N <- 1000
    for (strand in c("original", "reverse", "both")) {
        randomized <- vapply(rep(len, N), GENERATE_RANDOM_SEQ, FUN.VALUE="")
        barcodes <- sprintf(barcode.fmt, randomized)
        names(barcodes) <- seq_len(N)

        B <- DNAStringSet(barcodes)
        strandFUN <- CHOOSE_STRAND_FUN(strand)
        tmp <- tempfile(fileext=".fastq")
        writeXStringSet(strandFUN(B), filepath=tmp, format="fastq")

        out <- countRandomBarcodes(tmp, template=template, strand=strand)
        tab <- table(randomized)
        o <- order(names(tab))
        expect_identical(names(tab)[o], out$sequences)
        expect_identical(as.integer(unname(tab))[o], out$counts)

        # Same results when you stick a bunch of random crap to the start and end.
        STICKER(barcodes, tmp, out, strand=strand, strandFUN=strandFUN)
    }
})

test_that("matrixOfRandomBarcodes works as expected with names", {
    N <- 20
    randomized <- vapply(rep(len, N), GENERATE_RANDOM_SEQ, FUN.VALUE="")
    barcodes <- sprintf(barcode.fmt, randomized)
    names(barcodes) <- seq_len(N)

    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

    mat <- matrixOfRandomBarcodes(tmp, template=template)
    expect_identical(colnames(mat), basename(tmp))
    expect_identical(rownames(mat), rowData(mat)$sequences)
})
