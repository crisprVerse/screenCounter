# This tests the dual counting capability.
# library(testthat); library(screenCounter); source("setup.R"); source("test-dual.R")

set.seed(100000)
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

test_that("dual counting works as expected for edits", {
    barcodes <- c(
        "ACGTGGGGGGGGGGACGT",
        "ACGTGGGGCGGGGGACGT",
        "ACGTGGGGGGGGGACGT",
        "ACGTGGGGGGGGGGGACGT"  
    )
    names(barcodes) <- seq_along(barcodes)

    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

    choices <- strrep(BASES, 10)
    output <- countDualBarcodes(c(tmp, tmp), choices=DataFrame(choices, choices), template=template1)
    expect_identical(sum(output$counts), 1L)
    ref <- countSingleBarcodes(tmp, choices=choices, template=template1, strand="original")
    expect_identical(ref$counts, output$counts)

    # Throwing in scalar specifications.
    output <- countDualBarcodes(c(tmp, tmp), choices=DataFrame(choices, choices), template=template1, insertion=1)
    expect_identical(sum(output$counts), 2L)
    ref <- countSingleBarcodes(tmp, choices=choices, template=template1, strand="original", insertion=1)
    expect_identical(ref$counts, output$counts)
    
    output <- countDualBarcodes(c(tmp, tmp), choices=DataFrame(choices, choices), template=template1, substitution=1)
    expect_identical(sum(output$counts), 2L)
    ref <- countSingleBarcodes(tmp, choices=choices, template=template1, strand="original", substitution=1)
    expect_identical(ref$counts, output$counts)

    output <- countDualBarcodes(c(tmp, tmp), choices=DataFrame(choices, choices), template=template1, deletion=1)
    expect_identical(sum(output$counts), 2L)
    ref <- countSingleBarcodes(tmp, choices=choices, template=template1, strand="original", deletion=1)
    expect_identical(ref$counts, output$counts)

    # Works with vectors.
    barcodes <- c(
        "ACGTGGGGCGGGGGACGT",
        "ACGTGGGGGGGGGGACGT",
        "ACGTGGGGGGGGGGGACGT",
        "ACGTGGGGGGGGGACGT"
    )
    names(barcodes) <- seq_along(barcodes)

    tmp2 <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp2, format="fastq")

    output <- countDualBarcodes(c(tmp, tmp2), choices=DataFrame(choices, choices), template=template1)
    expect_identical(sum(output$counts), 0L)

    output <- countDualBarcodes(c(tmp, tmp2), choices=DataFrame(choices, choices), template=template1, 
        substitution=c(0, 1), insertion=c(0, 1), deletion=c(0, 1))
    expect_identical(sum(output$counts), 1L)

    output <- countDualBarcodes(c(tmp, tmp2), choices=DataFrame(choices, choices), template=template1, 
        substitution=c(1, 0), insertion=c(1, 0), deletion=c(1, 0))
    expect_identical(sum(output$counts), 1L)

    output <- countDualBarcodes(c(tmp, tmp2), choices=DataFrame(choices, choices), template=template1, 
        substitution=c(1, 1), insertion=c(1, 1), deletion=c(1, 1))
    expect_identical(sum(output$counts), 4L)
})

##############################################
##############################################

combos <- expand.grid(POOL1, POOL2)
choices <- DataFrame(X=as.character(combos[,1]), Y=as.character(combos[,2]))

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
    m <- match(DataFrame(X=POOL1[i], Y=POOL2[j]), choices)
    expect_identical(output$counts, tabulate(m, nbins=nrow(choices))) # NA's are dismissed.

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
    expect_identical(output$counts, tabulate(match(DataFrame(X=POOL1[i], Y=POOL2[j]), choices), nbins=nrow(choices))) 

    # Checking diagnostics:
    collated <- metadata(output)
    expect_identical(collated$none, sum(i > nbarcodes1 & j > nbarcodes2))
    expect_identical(collated$barcode1.only, sum(i <= nbarcodes1 & j > nbarcodes2))
    expect_identical(collated$barcode2.only, sum(i > nbarcodes1 & j <= nbarcodes2))
    expect_identical(collated$invalid.pair, 0L)
    expect_equal(collated$none + collated$barcode1.only + collated$barcode2.only + sum(output$counts), N)
})

test_that("dual counting handles randomization correctly", {
    N <- 5000L

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

    # Checking diagnostics:
    collated <- metadata(output)
    expect_identical(collated$none, N)
    collated <- metadata(output2)
    expect_identical(collated$none, 0L)
    expect_identical(collated$original.orientation, N)
    expect_identical(collated$reverse.orientation, N)

    # Trying with a symmetric construct.
    i <- sample(nbarcodes1, N, replace=TRUE)
    barcodes <- sprintf(barcode.fmt1, POOL1[i])
    names(barcodes) <- seq_len(N)
    barcodes2 <- sprintf(barcode.fmt2, POOL1[i])
    names(barcodes2) <- seq_len(N)

    writeXStringSet(DNAStringSet(c(barcodes, barcodes2)), filepath=tmp, format="fastq")
    writeXStringSet(DNAStringSet(c(barcodes2, barcodes)), filepath=tmp2, format="fastq")

    choices <- DataFrame(X=POOL1, Y=POOL1)
    output <- countDualBarcodes(c(tmp, tmp2), choices=choices, template=c(template1, template2))
    expect_identical(as.data.frame(output[,1:2]), as.data.frame(choices))
    expect_identical(output$counts, tabulate(i, nbins=length(POOL1)))
})

test_that("dual counting reports invalid pairs correctly", {
    N <- 5000

    i <- sample(nbarcodes1, N, replace=TRUE)
    barcodes <- sprintf(barcode.fmt1, POOL1[i])
    names(barcodes) <- seq_len(N)

    j <- sample(nbarcodes2, N, replace=TRUE)
    barcodes2 <- sprintf(barcode.fmt2, POOL2[j])
    names(barcodes2) <- seq_len(N)

    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")
    tmp2 <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes2), filepath=tmp2, format="fastq")

    ref <- countDualBarcodes(c(tmp, tmp2), choices=choices, template=c(template1, template2))
    sub <- countDualBarcodes(c(tmp, tmp2), choices=choices, template=c(template1, template2), include.invalid=TRUE)
    expect_true(all(sub$valid))
    sub$valid <- NULL
    expect_identical(ref, sub)

    # Mocking up a situation with only a subset of combinations, but all barcodes present.
    N <- max(length(POOL1), length(POOL2))
    subchoices <- DataFrame(X=factor(rep(POOL1, length.out=N)), Y=factor(rep(POOL2, length.out=N)))
    sub <- countDualBarcodes(c(tmp, tmp2), choices=subchoices, template=c(template1, template2), include.invalid=TRUE)

    expect_true(all(sub$valid[1:nrow(subchoices)]))
    expect_false(any(sub$valid[(nrow(subchoices)+1):nrow(sub)]))
    m <- match(ref[,1:2], sub[,1:2])
    nzero <- ref$counts!=0
    expect_identical(ref$counts[nzero], sub$counts[m[nzero]])
    expect_identical(sum(nzero), nrow(sub))

    # Working up a situation involving randomization.
    writeXStringSet(DNAStringSet(c(barcodes, barcodes2)), filepath=tmp, format="fastq")
    writeXStringSet(DNAStringSet(c(barcodes2, barcodes)), filepath=tmp2, format="fastq")

    output2 <- countDualBarcodes(c(tmp, tmp2), choices=subchoices, template=c(template1, template2), include.invalid=TRUE, randomized=TRUE)
    expect_identical(as.data.frame(sub[,1:2]), as.data.frame(output2[,1:2]))
    expect_identical(sub$counts*2L, output2$counts)
    expect_identical(sub$valid, output2$valid)
})

##############################################
##############################################

SPAWN_MULTI_FILES <- function(Nchoices) {
    collected <- list()

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

    collected
}

test_that("matrix summarization works as expected", {
    Nchoices <- c(5000, 2000, 1000)
    collected <- SPAWN_MULTI_FILES(Nchoices)

    se <- matrixOfDualBarcodes(collected, choices=choices, template=c(template1, template2))
    expect_true(all(c("none", "barcode1.only", "barcode2.only", "invalid.pair") %in% colnames(colData(se))))
    expect_equivalent(colSums(assay(se)), Nchoices)

    # Same result after including invalid pairs.
    se2 <- matrixOfDualBarcodes(collected, choices=choices, template=c(template1, template2), include.invalid=TRUE)
    expect_true(length(rowData(se2)$valid) && all(rowData(se2)$valid))
    rowData(se2)$valid <- NULL

    o1 <- order(rowData(se))
    o2 <- order(rowData(se2))
    expect_identical(se[o1,], se2[o2,])

    # Respects row names.
    rownames(choices) <- paste0("PAIR_", seq_len(nrow(choices)))
    full <- matrixOfDualBarcodes(collected, choices=choices, template=c(template1, template2))
    expect_identical(rownames(full), rownames(choices))
})

test_that("matrix summarization works as expected with invalid pairs", {
    Nchoices <- c(5000, 2000, 1000)
    collected <- SPAWN_MULTI_FILES(Nchoices)

    rownames(choices) <- paste0("PAIR_", seq_len(nrow(choices)))
    full <- matrixOfDualBarcodes(collected, choices=choices, template=c(template1, template2))
    expect_identical(rownames(full), rownames(choices))

    # Correctly respects the subset.
    max.len <- max(length(POOL1), length(POOL2))
    subchoices <- DataFrame(X=rep(POOL1, length.out=max.len), Y=rep(POOL2, length.out=max.len))
    keep <- match(subchoices, choices)
    rownames(subchoices) <- rownames(choices)[keep]

    ref <- matrixOfDualBarcodes(collected, choices=subchoices, template=c(template1, template2))
    expect_identical(rowData(full)[keep,], rowData(ref))
    expect_identical(assay(full)[keep,], assay(ref))
    expect_identical(rownames(ref), rownames(subchoices))

    # Correctly reports valid pairs when we include invalids.
    se <- matrixOfDualBarcodes(collected, choices=subchoices, template=c(template1, template2), include.invalid=TRUE)
    expect_true(nrow(se) > nrow(ref))
    sub <- se[rowData(se)$valid,]
    rowData(sub)$valid <- NULL

    m <- match(rowData(sub), rowData(ref))
    expect_identical(ref[m,], sub)

    # Same results if we apply an ordered input. This checks that the conversion
    # of indices to sequences is done correctly for the invalid pairs.
    o <- do.call(order, c(as.list(subchoices), list(decreasing=TRUE)))
    subchoices2 <- subchoices[o,]
    keep <- match(subchoices, choices)

    ref <- matrixOfDualBarcodes(collected, choices=subchoices, template=c(template1, template2), include.invalid=TRUE)
    se2 <- matrixOfDualBarcodes(collected, choices=subchoices2, template=c(template1, template2), include.invalid=TRUE)
    expect_identical(nrow(ref), nrow(se2))

    m <- match(rowData(se2)[,1:2], rowData(ref)[,1:2])
    expect_false(any(is.na(m)))
    expect_identical(ref[m,], se2)
})
