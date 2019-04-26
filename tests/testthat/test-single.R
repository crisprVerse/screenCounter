# This tests countSingleBarcodes().
# library(testthat); library(gp.sa.screen); source("setup.R"); source("test-single.R")

library(Biostrings)
nbarcodes <- 20
POOL <- vapply(rep(10, nbarcodes), GENERATE_RANDOM_SEQ, FUN.VALUE="")

barcode.fmt <- "ACGT%sACGT"
template <- sprintf(barcode.fmt, strrep("N", nchar(POOL[1])))

STICKER <- function(barcodes, fname, out, ..., strandFUN=identity) {
    ADD_FLANKS(barcodes, fname, strandFUN=strandFUN)
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

        out <- countSingleBarcodes(tmp, POOL, template=template)
        tab <- tabulate(i, nbins=nbarcodes)
        expect_identical(tab, out)

        # Same results when you stick a bunch of random crap to the start and end.
        STICKER(barcodes, tmp, out, choices=POOL)
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
    out <- countSingleBarcodes(tmp, choices, template=template, sub=TRUE)
    expect_identical(out, c(0L, 0L, 2L, 0L))

    STICKER(barcodes, tmp, out, choices=choices, sub=TRUE)

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
    out <- countSingleBarcodes(tmp, choices, template=template, sub=TRUE)
    expect_identical(out, c(1L, 1L))

    STICKER(barcodes, tmp, out, choices=choices, sub=TRUE)
})

test_that("countSingleBarcodes works as expected with deletions", {
    barcodes <- c(
        "ACGTCCCCCCCCCCACGT",
        "ACGTCCCCCCCCCACGT",
        "ACGTCCCCCCCCACGT",  # not matched, two deletions.
        "ACGTCCCCCCCCCCAGT", # not matched, error in a constant region.
        "ACTCCCCCCCCCCACGT"  # not matched, error in a constant region.
    )
    names(barcodes) <- seq_along(barcodes)

    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

    choices <- strrep(BASES, 10)
    out <- countSingleBarcodes(tmp, choices, template=template, del=TRUE)
    expect_identical(out, c(0L, 2L, 0L, 0L))

    STICKER(barcodes, tmp, out, choices=choices, del=TRUE)

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
    out <- countSingleBarcodes(tmp, choices, template=template, del=TRUE)
    expect_identical(out, c(1L, 2L))

    STICKER(barcodes, tmp, out, choices=choices, del=TRUE)
})

test_that("countSingleBarcodes works as expected with strands", {
    for (strand in c("original", "reverse", "both")) {
        i <- sample(nbarcodes, 1000, replace=TRUE)
        barcodes <- sprintf(barcode.fmt, POOL[i])
        names(barcodes) <- seq_along(i)

        B <- DNAStringSet(barcodes)
        strandFUN <- CHOOSE_STRAND_FUN(strand)
        tmp <- tempfile(fileext=".fastq")
        writeXStringSet(strandFUN(B), filepath=tmp, format="fastq")

        out <- countSingleBarcodes(tmp, POOL, template=template, strand=strand)
        tab <- tabulate(i, nbins=nbarcodes)
        expect_identical(tab, out)

        # Same results when you stick a bunch of random crap to the start and end.
        STICKER(barcodes, tmp, out, choices=POOL, strand=strand, strandFUN=strandFUN)
    }
})
