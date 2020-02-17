# This tests countComboBarcodes().
# library(testthat); library(screenCounter); source("setup.R"); source("test-combo-count.R")

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

        tab <- table(factor(i1, seq_along(POOL1)), factor(i2, seq_along(POOL2)))
        ref <- list(as.vector(row(tab)), as.vector(col(tab)), as.vector(tab)) 
        keep <- ref[[3]]!=0L
        ref <- lapply(ref, "[", keep)
        o <- do.call(order, ref)
        ref <- lapply(ref, "[", o)

        CHECK_OUTPUT(out, ref[[1]], ref[[2]], ref[[3]], N)

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
    out <- countComboBarcodes(tmp, template, choices, sub=TRUE, indices=TRUE)
    CHECK_OUTPUT(out, ref1=2L, ref2=3L, count=3L, N=length(barcodes))

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
    out <- countComboBarcodes(tmp, template, choices=choices, sub=TRUE, indices=TRUE)
    CHECK_OUTPUT(out, ref1=rep(1:2, each=2), ref2=rep(1:2, 2), count=rep(1L, 4), N=length(barcodes))

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
    out <- countComboBarcodes(tmp, template, choices, del=TRUE, indices=TRUE)
    CHECK_OUTPUT(out, ref1=2L, ref2=3L, count=3L, N=length(barcodes))

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
    out <- countComboBarcodes(tmp, template, choices=choices, del=TRUE, indices=TRUE)
    CHECK_OUTPUT(out, ref1=c(1L, 1L, 2L), ref2=c(1L, 2L, 1L), count=rep(1L, 3), N=length(barcodes))

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

        out <- countComboBarcodes(tmp, template, list(POOL1, POOL2), strand=strand, indices=TRUE)

        tab <- table(factor(i1, seq_along(POOL1)), factor(i2, seq_along(POOL2)))
        ref <- list(as.vector(row(tab)), as.vector(col(tab)), as.vector(tab)) 
        keep <- ref[[3]]!=0L
        ref <- lapply(ref, "[", keep)
        o <- do.call(order, ref)
        ref <- lapply(ref, "[", o)

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
    expect_identical(colnames(rowData(se)), c("X1", "X2"))
}) 
