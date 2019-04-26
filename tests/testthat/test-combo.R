# This tests countComboBarcodes().
# library(testthat); library(gp.sa.screen); source("setup.R"); source("test-combo.R")

library(Biostrings)
POOL1 <- vapply(rep(10, 20), GENERATE_RANDOM_SEQ, FUN.VALUE="")
POOL2 <- vapply(rep(8, 10), GENERATE_RANDOM_SEQ, FUN.VALUE="") 

barcode.fmt <- "ACGT%sACGT%sACGT"
template <- sprintf(barcode.fmt, strrep("N", nchar(POOL1[1])), strrep("N", nchar(POOL2[1])))

STICKER <- function(barcodes, fname, out, ..., strandFUN=identity) {
    ADD_FLANKS(barcodes, fname, strandFUN=strandFUN)
    out2 <- countComboBarcodes(fname, template, ...)
    expect_identical(out, out2)
}

test_that("countComboBarcodes works as expected in basic mode", {
    for (N in c(1, 10, 100, 1000)) {
        i1 <- sample(length(POOL1), N, replace=TRUE)
        i2 <- sample(length(POOL2), N, replace=TRUE)
        barcodes <- sprintf(barcode.fmt, POOL1[i1], POOL2[i2])
        names(barcodes) <- seq_len(N)

        tmp <- tempfile(fileext=".fastq")
        writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

        out <- countComboBarcodes(tmp, template, list(POOL1, POOL2))

        tab <- table(factor(i1, seq_along(POOL1)), factor(i2, seq_along(POOL2)))
        ref <- list(as.vector(row(tab)), as.vector(col(tab)), as.vector(tab)) 
        keep <- ref[[3]]!=0L
        ref <- lapply(ref, "[", keep)
        o <- do.call(order, ref)
        ref <- lapply(ref, "[", o)

        expect_identical(out$combination[,1], ref[[1]])
        expect_identical(out$combination[,2], ref[[2]])
        expect_identical(out$count, ref[[3]])

        # Same results when you stick a bunch of random crap to the start and end.
        STICKER(barcodes, tmp, out, choices=list(POOL1, POOL2))
    }
})

test_that("countComboBarcodes works as expected with substitutions", {
    barcodes <- c(
        "ACGTCCCCCCCCCCACGTGGGGGGGGACGT",
        "ACGTCCCCCGCCCCACGTGGGGGGGGACGT",
        "ACGTCCCCCCCCCCACGTGGGGCGGGACGT",
        "ACGTCCCCGCCCCCACGTGGGGCGGGACGT", # not matched, error in both variable regions.
        "ACGTCCCCCCCCCCACGTGGGGGGGGAGGT", # not matched, error in a constant region.
        "ACGTCCCCCCCCCCATGTGGGGGGGGACGT", # not matched, error in a constant region.
        "CCGTCCCCCCCCCCACGTGGGGGGGGACGT"  # not matched, error in a constant region.
    )
    names(barcodes) <- seq_along(barcodes)

    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

    choices <- list(strrep(BASES, 10), strrep(BASES, 8))
    out <- countComboBarcodes(tmp, template, choices, sub=TRUE)
    expect_identical(out$combination[,1], 2L)
    expect_identical(out$combination[,2], 3L)
    expect_identical(out$count, 3L)

    STICKER(barcodes, tmp, out, choices=choices, sub=TRUE)

    # Handles conflicts correctly.
    barcodes <- c(
        "ACGTCCCCCCCCCCACGTGGGGGGGGACGT", 
        "ACGTCCCCCCCCCCACGTTGGGGGGGACGT",
        "ACGTCCCCCCCCCAACGTGGGGGGGGACGT",
        "ACGTCCCCCCCCCAACGTTGGGGGGGACGT",
        "ACGTCCCCCCCCCGACGTGGGGGGGGACGT", # ambiguous and removed.
        "ACGTCCCCCCCCCTACGTGGGGGGGGACGT", # ambiguous and removed.
        "ACGTCCCCCCCCCCACGTAGGGGGGGACGT", # ambiguous and removed.
        "ACGTCCCCCCCCCCACGTCGGGGGGGACGT"  # ambiguous and removed. 
    )
    names(barcodes) <- seq_along(barcodes)

    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

    choices <- list(c("CCCCCCCCCC", "CCCCCCCCCA"), c("GGGGGGGG", "TGGGGGGG"))
    out <- countComboBarcodes(tmp, template, choices=choices, sub=TRUE)
    expect_identical(out$combination[,1], rep(1:2, each=2))
    expect_identical(out$combination[,2], rep(1:2, 2))
    expect_identical(out$count, rep(1L, 4))

    STICKER(barcodes, tmp, out, choices=choices, sub=TRUE)
})

test_that("countComboBarcodes works as expected with deletions", {
    barcodes <- c(
        "ACGTCCCCCCCCCCACGTGGGGGGGGACGT",
        "ACGTCCCCCCCCCACGTGGGGGGGGACGT",
        "ACGTCCCCCCCCCCACGTGGGGGGGACGT",
        "ACGTCCCCCCCCCACGTGGGGGGGACGT",  # not matched, error in both variable regions.
        "ACGTCCCCCCCCCCACGTGGGGGGGGAGT", # not matched, error in a constant region.
        "ACGTCCCCCCCCCCAGTGGGGGGGGACGT", # not matched, error in a constant region.
        "ACGCCCCCCCCCCACGTGGGGGGGGACGT"  # not matched, error in a constant region.
    )
    names(barcodes) <- seq_along(barcodes)

    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

    choices <- list(strrep(BASES, 10), strrep(BASES, 8))
    out <- countComboBarcodes(tmp, template, choices, del=TRUE)
    expect_identical(out$combination[,1], 2L)
    expect_identical(out$combination[,2], 3L)
    expect_identical(out$count, 3L)

    STICKER(barcodes, tmp, out, choices=choices, del=TRUE)

    # Handles conflicts correctly.
    barcodes <- c(
        "ACGTCCCCCCCCCCACGTGGGGGGGGACGT", 
        "ACGTCCCCCCCCAACGTGGGGGGGGACGT", 
        "ACGTCCCCCCCCCCACGTTGGGGGGACGT", 
        "ACGTCCCCCCCCCACGTGGGGGGGGACGT", # ambiguous and removed.
        "ACGTCCCCCCCCCCACGTGGGGGGGACGT"  # ambiguous and removed. 
    )
    names(barcodes) <- seq_along(barcodes)

    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

    choices <- list(c("CCCCCCCCCC", "CCCCCCCCCA"), c("GGGGGGGG", "TGGGGGGG"))
    out <- countComboBarcodes(tmp, template, choices=choices, del=TRUE)
    expect_identical(out$combination[,1], c(1L, 1L, 2L))
    expect_identical(out$combination[,2], c(1L, 2L, 1L))
    expect_identical(out$count, rep(1L, 3))

    STICKER(barcodes, tmp, out, choices=choices, del=TRUE)
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

        out <- countComboBarcodes(tmp, template, list(POOL1, POOL2), strand=strand)

        tab <- table(factor(i1, seq_along(POOL1)), factor(i2, seq_along(POOL2)))
        ref <- list(as.vector(row(tab)), as.vector(col(tab)), as.vector(tab)) 
        keep <- ref[[3]]!=0L
        ref <- lapply(ref, "[", keep)
        o <- do.call(order, ref)
        ref <- lapply(ref, "[", o)

        expect_identical(out$combination[,1], ref[[1]])
        expect_identical(out$combination[,2], ref[[2]])
        expect_identical(out$count, ref[[3]])

        # Same results when you stick a bunch of random crap to the start and end.
        STICKER(barcodes, tmp, out, choices=list(POOL1, POOL2), strand=strand, strandFUN=strandFUN)
    }
})

