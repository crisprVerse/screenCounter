# This tests the dual counting capability.
# library(testthat); library(screenCounter); source("setup.R"); source("test-dual.R")

library(Biostrings)
nbarcodes1 <- 20
POOL1 <- vapply(rep(10, nbarcodes1), GENERATE_RANDOM_SEQ, FUN.VALUE="")
nbarcodes2 <- 50
POOL2 <- vapply(rep(15, nbarcodes2), GENERATE_RANDOM_SEQ, FUN.VALUE="")

barcode.fmt1 <- "ACGT%sACGT"
template1 <- sprintf(barcode.fmt1, strrep("N", nchar(POOL1[1])))
barcode.fmt2 <- "AGGA%sAGGA"
template2 <- sprintf(barcode.fmt2, strrep("N", nchar(POOL2[1])))

test_that("dual counting gives the same results as the single counter", {
    N <- 200

    # Vanilla approach with the same template.
    i <- sample(nbarcodes1, N, replace=TRUE)
    barcodes <- sprintf(barcode.fmt1, POOL1[i])
    names(barcodes) <- seq_len(N)

    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

    output <- countDualBarcodes(c(tmp, tmp), choices=DataFrame(POOL1, POOL1), template=template1)
    ref <- countSingleBarcodes(tmp, choices=POOL1, template=template1, strand="original")
    expect_identical(output$counts, ref$counts)

    # Checking that it responds to different templates and strands.
    tmp2 <- tempfile(fileext=".fastq")
    barcodes2 <- sprintf(barcode.fmt2, POOL1[i])
    names(barcodes2) <- seq_len(N)
    writeXStringSet(reverseComplement(DNAStringSet(barcodes2)), filepath=tmp2, format="fastq")

    output2 <- countDualBarcodes(c(tmp, tmp2), choices=DataFrame(POOL1, POOL1), 
        template=c(template1, template2), strand=c("original", "reverse"))
    expect_identical(output2$counts, ref$counts)

    output.x <- countDualBarcodes(c(tmp, tmp2), choices=DataFrame(POOL1, POOL1), template=template1)
    expect_false(identical(output.x$counts, ref$counts))
})

combos <- expand.grid(POOL1, POOL2)
choices <- DataFrame(X=combos[,1], Y=combos[,2])

test_that("dual counting works as expected for combinations", {
    N <- 5000

    # Vanilla example works as expected.
    i <- sample(nbarcodes1, N, replace=TRUE)
    barcodes <- sprintf(barcode.fmt1, POOL1[i])
    names(barcodes) <- seq_len(N)
    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

    j <- sample(nbarcodes2, N, replace=TRUE)
    barcodes2 <- sprintf(barcode.fmt2, POOL2[j])
    names(barcodes2) <- seq_len(N)
    tmp2 <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes2), filepath=tmp2, format="fastq")

    output <- countDualBarcodes(c(tmp, tmp2), choices=choices, template=c(template1, template2))
    expect_identical(as.data.frame(output[,1:2]), as.data.frame(choices))
    expect_identical(output$counts, tabulate(match(DataFrame(X=POOL1[i], Y=POOL2[j]), choices)))

    # Works when only a subset of the combinations are valid.
    keep <- sample(nrow(choices), nrow(choices)/2)
    choices2 <- choices[keep,]
    output2 <- countDualBarcodes(c(tmp, tmp2), choices=choices2, template=c(template1, template2))
    expect_identical(as.data.frame(output2[,1:2]), as.data.frame(choices[keep,]))
    expect_identical(output2$counts, output$counts[keep])
    expect_identical(metadata(output2)$invalid.pair, as.integer(N) - sum(output2$counts))
})

test_that("dual counting reports diagnostic values correctly", {
    N <- 5000

    # Adding about 50% missing values.
    POOL1.0 <- c(POOL1, rep(strrep("A", nchar(POOL1[1])), nbarcodes1))
    POOL2.0 <- c(POOL2, rep(strrep("A", nchar(POOL2[1])), nbarcodes2))

    # Vanilla example works as expected.
    i <- sample(length(POOL1.0), N, replace=TRUE)
    barcodes <- sprintf(barcode.fmt1, POOL1.0[i])
    names(barcodes) <- seq_len(N)
    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

    j <- sample(length(POOL2.0), N, replace=TRUE) 
    barcodes2 <- sprintf(barcode.fmt2, POOL2.0[j])
    names(barcodes2) <- seq_len(N)
    tmp2 <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes2), filepath=tmp2, format="fastq")

    output <- countDualBarcodes(c(tmp, tmp2), choices=choices, template=c(template1, template2))
    expect_identical(as.data.frame(output[,1:2]), as.data.frame(choices))
    expect_identical(output$counts, tabulate(match(DataFrame(X=POOL1[i], Y=POOL2[j]), choices))) # NA's are dismissed.

    # Checking diagnostics:
    collated <- metadata(output)
    expect_identical(collated$none, sum(i > nbarcodes1 & j > nbarcodes2))
    expect_identical(collated$barcode1.only, sum(i <= nbarcodes1 & j > nbarcodes2))
    expect_identical(collated$barcode2.only, sum(i > nbarcodes1 & j <= nbarcodes2))
    expect_identical(collated$invalid.pair, 0L)
    expect_equal(collated$none + collated$barcode1.only + collated$barcode2.only + sum(output$counts), N)
})

test_that("dual counting handles randomization correctly", {
    N <- 5000

    i <- sample(nbarcodes1, N, replace=TRUE)
    barcodes <- sprintf(barcode.fmt1, POOL1[i])
    names(barcodes) <- seq_len(N)

    j <- sample(nbarcodes2, N, replace=TRUE)
    barcodes2 <- sprintf(barcode.fmt2, POOL2[j])
    names(barcodes2) <- seq_len(N)

    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(c(barcodes, barcodes2)), filepath=tmp, format="fastq")
    tmp2 <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(c(barcodes2, barcodes)), filepath=tmp2, format="fastq")

    output <- countDualBarcodes(c(tmp, tmp2), choices=choices, template=c(template1, template2))
    output2 <- countDualBarcodes(c(tmp, tmp2), choices=choices, template=c(template1, template2), randomized=TRUE)

    expect_identical(as.data.frame(output[,1:2]), as.data.frame(output2[,1:2]))
    expect_identical(output$counts*2L, output2$counts)
})

test_that("matrix summarization works as expected", {
    collected <- list()
    Nchoices <- c(5000, 2000, 1000)

    for (x in seq_along(Nchoices)) {
        N <- Nchoices[x]

        # Vanilla example works as expected.
        i <- sample(nbarcodes1, N, replace=TRUE)
        barcodes <- sprintf(barcode.fmt1, POOL1[i])
        names(barcodes) <- seq_len(N)
        tmp <- tempfile(fileext=".fastq")
        writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

        j <- sample(nbarcodes2, N, replace=TRUE)
        barcodes2 <- sprintf(barcode.fmt2, POOL2[j])
        names(barcodes2) <- seq_len(N)
        tmp2 <- tempfile(fileext=".fastq")
        writeXStringSet(DNAStringSet(barcodes2), filepath=tmp2, format="fastq")

        collected[[x]] <- c(tmp, tmp2)
    }

    se <- matrixOfDualBarcodes(collected, choices=choices, template=c(template1, template2))
    expect_true(all(c("none", "barcode1.only", "barcode2.only", "invalid.pair") %in% colnames(colData(se))))
    expect_equivalent(colSums(assay(se)), Nchoices)
})
