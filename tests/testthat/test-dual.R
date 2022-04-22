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

    template2b <- sprintf(barcode.fmt2, strrep("N", nchar(POOL1[1])))
    output2 <- countDualBarcodes(c(tmp, tmp2), choices=DataFrame(POOL1, POOL1), 
        template=c(template1, template2b), strand=c("original", "reverse"))
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
    output <- countDualBarcodes(c(tmp, tmp), choices=DataFrame(choices, choices), template=template1, substitution=1)
    expect_identical(sum(output$counts), 2L)
    ref <- countSingleBarcodes(tmp, choices=choices, template=template1, strand="original", substitution=1)
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

    output <- countDualBarcodes(c(tmp, tmp2), choices=DataFrame(choices, choices), template=template1, substitution=c(0, 1)) 
    expect_identical(sum(output$counts), 1L)

    output <- countDualBarcodes(c(tmp, tmp2), choices=DataFrame(choices, choices), template=template1, substitution=c(1, 0))
    expect_identical(sum(output$counts), 1L)

    output <- countDualBarcodes(c(tmp, tmp2), choices=DataFrame(choices, choices), template=template1, substitution=c(1, 1))
    expect_identical(sum(output$counts), 2L)
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
    m <- S4Vectors::match(DataFrame(X=POOL1[i], Y=POOL2[j]), choices)
    expect_identical(output$counts, tabulate(m, nbins=nrow(choices))) # NA's are dismissed.

    # Works when only a subset of the combinations are valid.
    keep <- sample(nrow(choices), nrow(choices)/2)
    choices2 <- choices[keep,]
    output2 <- countDualBarcodes(c(tmp, tmp2), choices=choices2, template=c(template1, template2))
    expect_identical(as.data.frame(output2[,1:2]), as.data.frame(choices[keep,]))
    expect_identical(output2$counts, output$counts[keep])
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

    # Trying with a symmetric construct.
    i <- sample(nbarcodes1, N, replace=TRUE)
    barcodes <- sprintf(barcode.fmt1, POOL1[i])
    names(barcodes) <- seq_len(N)
    barcodes2 <- sprintf(barcode.fmt2, POOL1[i])
    names(barcodes2) <- seq_len(N)

    writeXStringSet(DNAStringSet(c(barcodes, barcodes2)), filepath=tmp, format="fastq")
    writeXStringSet(DNAStringSet(c(barcodes2, barcodes)), filepath=tmp2, format="fastq")

    choices <- DataFrame(X=POOL1, Y=POOL1)
    template2b <- sprintf(barcode.fmt2, strrep("N", nchar(POOL1[1])))
    output <- countDualBarcodes(c(tmp, tmp2), choices=choices, template=c(template1, template2b))
    expect_identical(as.data.frame(output[,1:2]), as.data.frame(choices))
    expect_identical(output$counts, tabulate(i, nbins=length(POOL1)))
})

test_that("dual counting handles randomization edge cases", {
    barcodes <- c(
        "AAAAAAAAA",
        "AAAAAAAAA",
        "AAAAAACAA",
        "AAAAAACAA"
    )
    names(barcodes) <- seq_along(barcodes)
    tmp <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")

    barcodes2 <- c(
        "AAAAAAAAA",
        "AAAAAACAA",
        "AAAAAAAAA",
        "AAAAAACAA"
    )
    names(barcodes2) <- seq_along(barcodes2)
    tmp2 <- tempfile(fileext=".fastq")
    writeXStringSet(DNAStringSet(barcodes2), filepath=tmp2, format="fastq")

    template <- "---------"
    output <- countDualBarcodes(c(tmp, tmp2), choices=DataFrame(first="AAAAAAAAA", second="AAAAAAAAA"), template=template)
    expect_identical(output$counts, 1L)

    output2 <- countDualBarcodes(c(tmp, tmp2), choices=DataFrame(first="AAAAAAAAA", second="AAAAAAAAA"), template=template, randomized=TRUE)
    expect_identical(as.data.frame(output), as.data.frame(output2))

    output3 <- countDualBarcodes(c(tmp, tmp2), choices=DataFrame(first="AAAAAAAAA", second="AAAAAAAAA"), template=template, substitutions=c(1, 0))
    expect_identical(output3$counts, 2L)

    output4 <- countDualBarcodes(c(tmp, tmp2), choices=DataFrame(first="AAAAAAAAA", second="AAAAAAAAA"), template=template, substitutions=c(0, 1))
    expect_identical(output4$counts, 2L)

    output5 <- countDualBarcodes(c(tmp, tmp2), choices=DataFrame(first="AAAAAAAAA", second="AAAAAAAAA"), template=template, substitutions=c(0, 1), randomized=TRUE)
    expect_identical(output5$counts, 3L) # doesn't match the case with two subs 

    output6a <- countDualBarcodes(c(tmp, tmp2), choices=DataFrame(first="AAAAAAAAA", second="AAAAAAAAA"), template=template, substitutions=1)
    output6b <- countDualBarcodes(c(tmp, tmp2), choices=DataFrame(first="AAAAAAAAA", second="AAAAAAAAA"), template=template, substitutions=1, randomized=TRUE)
    expect_identical(as.data.frame(output6a), as.data.frame(output6b))

    # Competition between matches in both orientations is resolved correctly.
    output7 <- countDualBarcodes(c(tmp, tmp2), choices=DataFrame(first=c("AAAAAAAAA", "AAAAAACAA"), second="AAAAAAAAA"), template=template, substitutions=1, randomized=TRUE)
    expect_identical(output7$counts[1], 2L) # 1 and 2 (mismatch is accepted by the use-first policy).
    expect_identical(output7$counts[2], 2L) # 3 and 4.
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
    metadata(sub) <- metadata(ref)
    expect_identical(ref, sub)

    # Mocking up a situation with only a subset of combinations, but all barcodes present.
    N <- max(length(POOL1), length(POOL2))
    subchoices <- DataFrame(X=factor(rep(POOL1, length.out=N)), Y=factor(rep(POOL2, length.out=N)))
    sub <- countDualBarcodes(c(tmp, tmp2), choices=subchoices, template=c(template1, template2), include.invalid=TRUE)

    expect_true(all(sub$valid[1:nrow(subchoices)]))
    expect_false(any(sub$valid[(nrow(subchoices)+1):nrow(sub)]))
    m <- S4Vectors::match(ref[,1:2], sub[,1:2])
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
    expect_true(all(c("npairs") %in% colnames(colData(se))))
    expect_equivalent(colSums(assay(se)), Nchoices)

    # Same result after including invalid pairs.
    se2 <- matrixOfDualBarcodes(collected, choices=choices, template=c(template1, template2), include.invalid=TRUE)
    expect_true(length(rowData(se2)$valid) && all(rowData(se2)$valid))
    rowData(se2)$valid <- NULL

    o1 <- order(rowData(se))
    o2 <- order(rowData(se2))
    colData(se2) <- colData(se2)[,c("paths1", "paths2", "npairs")]
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
    keep <- S4Vectors::match(subchoices, choices)
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

    m <- S4Vectors::match(rowData(sub), rowData(ref))
    colData(sub) <- colData(sub)[,c("paths1", "paths2", "npairs")]
    expect_identical(ref[m,], sub)

    # Same results if we apply an ordered input. This checks that the conversion
    # of indices to sequences is done correctly for the invalid pairs.
    o <- do.call(order, c(as.list(subchoices), list(decreasing=TRUE)))
    subchoices2 <- subchoices[o,]
    keep <- S4Vectors::match(subchoices, choices)

    ref <- matrixOfDualBarcodes(collected, choices=subchoices, template=c(template1, template2), include.invalid=TRUE)
    se2 <- matrixOfDualBarcodes(collected, choices=subchoices2, template=c(template1, template2), include.invalid=TRUE)
    expect_identical(nrow(ref), nrow(se2))

    m <- S4Vectors::match(rowData(se2)[,1:2], rowData(ref)[,1:2])
    expect_false(any(is.na(m)))
    expect_identical(ref[m,], se2)
})
