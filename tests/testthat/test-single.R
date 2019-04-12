# This tests countSingleBarcodes().
# library(testthat); library(gp.sa.screen); source("setup.R"); source("test-single.R")

library(Biostrings)
nbarcodes <- 20
POOL <- vapply(rep(10, nbarcodes), GENERATE_RANDOM_SEQ, FUN.VALUE="")

barcode.fmt <- "ACGT%sACGT"
template <- sprintf(barcode.fmt, strrep("N", nchar(POOL[1])))

STICKER <- function(barcodes, fname, out, ...) {
    ADD_FLANKS(barcodes, fname)
    out2 <- countSingleBarcodes(fname, template=template, ...)
    expect_identical(out, out2)
}

test_that("countSingleBarcodes works as expected in basic mode", {
    for (N in c(1, 10, 100, 1000)) {
        i <- sample(nbarcodes, N, replace=TRUE)
        barcodes <- sprintf(barcode.fmt, POOL[i])
        names(barcodes) <- seq_len(N)

        tmp <- tempfile(fileext=".fastq")
        writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

        out <- countSingleBarcodes(tmp, POOL, template=template, sub=FALSE, del=FALSE)
        tab <- tabulate(i, nbins=nbarcodes)
        expect_identical(tab, out)

        # Same results when you stick a bunch of random crap to the start and end.
        STICKER(barcodes, tmp, out, sub=FALSE, choices=POOL, del=FALSE)
    }
})

test_that("countSingleBarcodes works as expected with substitutions", {
    barcodes <- c(
        "ACGTGGGGGGGGGGACGT",
        "ACGTGGGGCGGGGGACGT",
        "ACGTGGGGCCGGGGACGT", # not matched, two substitutions.
        "ACGTGGGGGGGGGGATGT", # not matched, error in a constant region.
        "CCGTGGGGGGGGGGACGT"  # not matched, error in a constant region.
    )
    names(barcodes) <- seq_along(barcodes)

    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

    choices <- strrep(BASES, 10)
    out <- countSingleBarcodes(tmp, choices, template=template, del=FALSE)
    expect_identical(out, c(0L, 0L, 2L, 0L))

    STICKER(barcodes, tmp, out, choices=choices, del=FALSE)

    # Handles conflicts correctly.
    barcodes <- c(
        "ACGTCCCCCCCCCCACGT", 
        "ACGTCCCCCCCCCAACGT",
        "ACGTCCCCCCCCCGACGT", # ambiguous and removed.
        "ACGTCCCCCCCCCTACGT"  # ambiguous and removed.
    )
    names(barcodes) <- seq_along(barcodes)

    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

    choices <- c("CCCCCCCCCC", "CCCCCCCCCA")
    out <- countSingleBarcodes(tmp, choices, template=template, del=FALSE)
    expect_identical(out, c(1L, 1L))

    STICKER(barcodes, tmp, out, choices=choices, del=FALSE)
})

test_that("countSingleBarcodes works as expected with deletions", {
    barcodes <- c(
        "ACGTCCCCCCCCCCACGT",
        "ACGTCCCCCCCCCACGT",
        "ACGTCCCCCCCCACGT",  # not matched, two deletions.
        "ACGTCCCCCCCCCCAGT", # not matched, error in a constant region.
        "CGTCCCCCCCCCCACGT"  # not matched, error in a constant region.
    )
    names(barcodes) <- seq_along(barcodes)

    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

    choices <- strrep(BASES, 10)
    out <- countSingleBarcodes(tmp, choices, template=template, sub=FALSE)
    expect_identical(out, c(0L, 2L, 0L, 0L))

    STICKER(barcodes, tmp, out, choices=choices, sub=FALSE)

    # Handles conflicts correctly.
    barcodes <- c(
        "ACGTCCCCCCCCCCACGT", 
        "ACGTCCCCCCCCCAACGT", 
        "ACGTCCCCCCCCAACGT", 
        "ACGTCCCCCCCCCACGT"  # ambiguous and removed.
    )
    names(barcodes) <- seq_along(barcodes)

    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

    choices <- c("CCCCCCCCCC", "CCCCCCCCCA")
    out <- countSingleBarcodes(tmp, choices, template=template, sub=FALSE)
    expect_identical(out, c(1L, 2L))

    STICKER(barcodes, tmp, out, choices=choices, sub=FALSE)
})


