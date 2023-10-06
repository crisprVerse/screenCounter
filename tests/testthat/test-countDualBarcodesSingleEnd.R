# This tests the dual counting capability.
# library(testthat); library(screenCounter); source("setup.R"); source("test-countDualBarcodesSingleEnd.R")

set.seed(100000)
library(Biostrings)
nbarcodes <- 20
POOL1 <- vapply(rep(10, nbarcodes), GENERATE_RANDOM_SEQ, FUN.VALUE="")
POOL2 <- vapply(rep(15, nbarcodes), GENERATE_RANDOM_SEQ, FUN.VALUE="")

base.format <- "ACGT%sTGCAAGGA%sAGGA"
template <- sprintf(base.format, strrep("N", nchar(POOL1[1])), strrep("N", nchar(POOL2[1])))

test_that("countDualBarcodesSingleEnd gives the same results as the combinatorial counter", {
    N <- 200

    # Vanilla approach with the same template.
    i <- sample(nbarcodes, N, replace=TRUE)
    barcodes <- sprintf(base.format, POOL1[i], POOL2[i])
    names(barcodes) <- seq_len(N)

    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

    output <- countDualBarcodesSingleEnd(tmp, template=template, choices=DataFrame(POOL1, POOL2))
    ref <- countComboBarcodes(tmp, template=template, choices=List(POOL1, POOL2), strand="original")
    expect_identical(output$counts, ref$counts)

    # Checking that it responds to a different strand.
    tmp2 <- tempfile(fileext=".fastq")
    writeXStringSet(reverseComplement(DNAStringSet(barcodes)), filepath=tmp2, format="fastq")

    output2 <- countDualBarcodesSingleEnd(tmp2, template=template, choices=DataFrame(POOL1, POOL2))
    expect_identical(output2$counts, ref$counts)
    output2 <- countDualBarcodesSingleEnd(tmp2, template=template, choices=DataFrame(POOL1, POOL2), strand="original")
    expect_identical(sum(output2$counts), 0L)
})

test_that("countDualBarcodesSingleEnd works as expected for edits", {
    barcodes <- c(
        "ACGTGGGGGGGGGGTGCAAGGAAAAAAAAAAAAAAAAAGGA",
        "ACGTGGGGGGGGGGTGCAAGGAAAAAAAAAAATAAAAAGGA",
        "ACGTGGGGCGGGGGTGCAAGGAAAAAAAAAAATAAAAAGGA"
    )
    names(barcodes) <- seq_along(barcodes)

    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

    choices1 <- strrep(BASES, 10)
    choices2 <- strrep("A", 15)
    output <- countDualBarcodesSingleEnd(tmp, choices=DataFrame(choices1, choices2), template=template)
    expect_identical(output$counts, c(0L, 0L, 1L, 0L))

    output <- countDualBarcodesSingleEnd(tmp, choices=DataFrame(choices1, choices2), template=template, substitutions=1)
    expect_identical(output$counts, c(0L, 0L, 2L, 0L))

    output <- countDualBarcodesSingleEnd(tmp, choices=DataFrame(choices1, choices2), template=template, substitutions=2)
    expect_identical(output$counts, c(0L, 0L, 3L, 0L))
})

test_that("countDualBarcodesSingleEnd works when finding the best match", {
    N <- 200

    # Adding some mutations.
    i <- sample(nbarcodes, N, replace=TRUE)
    barcodes <- sprintf(base.format, POOL1[i], POOL2[i])
    for (j in seq_along(barcodes)) {
        current <- barcodes[j]
        end <- nchar(current)
        pos <- sample(end, 1)
        barcodes[j] <- paste0(substr(current, 1, pos - 1L), "N", substr(current, pos + 1L, end))
    }

    i2 <- sample(nbarcodes, N, replace=TRUE)
    barcodes2 <- sprintf(base.format, POOL1[i2], POOL2[i2])
    final.barcodes <- paste0(barcodes, "acgatcgatcgatcga", barcodes2)
    names(final.barcodes) <- seq_len(N)

    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(final.barcodes), filepath=tmp, format="fastq")

    ref <- countDualBarcodesSingleEnd(tmp, choices=DataFrame(POOL1, POOL2), template=template, substitutions=1, find.best=FALSE)
    expect_identical(ref$counts, tabulate(i, nbarcodes))

    output <- countDualBarcodesSingleEnd(tmp, choices=DataFrame(POOL1, POOL2), template=template, substitutions=1, find.best=TRUE)
    expect_identical(output$counts, tabulate(i2, nbarcodes))
})

test_that("dual counting reports invalid pairs correctly", {
    N <- 5000

    i1 <- sample(nbarcodes, N, replace=TRUE)
    i2 <- sample(nbarcodes, N, replace=TRUE)
    barcodes <- sprintf(base.format, POOL1[i1], POOL2[i2])
    names(barcodes) <- seq_len(N)

    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

    raw <- countDualBarcodesSingleEnd(tmp, choices=DataFrame(POOL1, POOL2), template=template)
    keep <- i1 == i2
    expect_identical(raw$counts, tabulate(i1[keep], nbarcodes))

    output <- countDualBarcodesSingleEnd(tmp, choices=DataFrame(POOL1, POOL2), template=template, include.invalid=TRUE)
    expect_identical(as.data.frame(raw), as.data.frame(output[output$valid,1:3]))

    everything <- countComboBarcodes(tmp, choices=List(POOL1, POOL2), template=template)
    m <- BiocGenerics::match(everything$combinations, DataFrame(first=output[,1], second=output[,2]))
    expect_identical(nrow(everything), nrow(output))
    expect_false(anyNA(m))
    expect_false(anyDuplicated(m) > 0)
    expect_identical(everything$counts, output$counts[m])
})

test_that("matrix summarization works as expected", {
    Nchoices <- c(500, 200, 100)
    fpaths <- character(0)
    selected <- list()

    for (N in Nchoices) {
        i1 <- sample(nbarcodes, N, replace=TRUE)
        i2 <- sample(nbarcodes, N, replace=TRUE)
        barcodes <- sprintf(base.format, POOL1[i1], POOL2[i2])
        names(barcodes) <- seq_len(N)

        tmp <- tempfile(fileext=".fastq")
        writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")
        fpaths <- c(fpaths, tmp)

        selected <- c(selected, list(list(i1, i2)))
    }

    # Ignoring invalid reads.
    se <- matrixOfDualBarcodesSingleEnd(fpaths, choices=DataFrame(POOL1, POOL2), template=template)
    expect_equal(se$nreads, Nchoices)

    counts <- assay(se)
    for (i in seq_len(ncol(counts))) {
        current <- selected[[i]]
        keep <- current[[1]] == current[[2]]
        expect_identical(counts[,i], tabulate(current[[1]][keep], nbarcodes))
    }

    # Including invalid reads.
    se2 <- matrixOfDualBarcodesSingleEnd(fpaths, choices=DataFrame(POOL1, POOL2), template=template, include.invalid=TRUE)
    o <- order(rowData(se))
    counts2 <- assay(se2)
    expect_identical(counts2[rowData(se2)$valid,], counts[o,])

    for (i in seq_len(ncol(counts))) {
        current <- selected[[i]]
        is.invalid <- current[[1]] != current[[2]]
        expect_identical(se2$invalid.reads[i], sum(is.invalid))
        expect_identical(sum(counts2[!rowData(se2)$valid,i]), sum(is.invalid))
    }
})
