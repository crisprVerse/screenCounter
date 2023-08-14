# This tests that the matchBarcodes function works as expected.
# library(testthat); library(screenCounter); source("test-matchBarcodes.R")

test_that("matchBarcodes works as expected in simple cases", {
    choices <- c("AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT")

    out <- matchBarcodes(c("AAAAAA", "AAATAA"), choices)
    expect_identical(out$index, c(1L, NA))
    expect_identical(out$mismatches, c(0L, NA))

    out <- matchBarcodes(c("AAAAAA", "AAATAA"), choices, substitutions=1)
    expect_identical(out$index, c(1L, 1L))
    expect_identical(out$mismatches, c(0L, 1L))

    out <- matchBarcodes(c("AAAAAA", "AAATAA"), choices, reverse=TRUE)
    expect_identical(out$index, c(4L, NA))
    expect_identical(out$mismatches, c(0L, NA))

    out <- matchBarcodes(c("AAAAAA", "AAATAA"), choices, reverse=TRUE, substitutions=1)
    expect_identical(out$index, c(4L, 4L))
    expect_identical(out$mismatches, c(0L, 1L))
})

test_that("matchBarcodes works as expected for IUPAC codes", {
    choices <- c("AAARAA", "CCCYCC", "GGGMGG", "TTTSTT")

    out <- matchBarcodes(c("AAAAAA", "AAAGAA"), choices)
    expect_identical(out$index, c(1L, 1L))
    expect_identical(out$mismatches, c(0L, 0L))

    out <- matchBarcodes(c("AAAAAA", "AAAGAA", "AAGAAA"), choices, reverse=TRUE)
    expect_identical(out$index, c(NA_integer_, NA_integer_, 4L))
    expect_identical(out$mismatches, c(NA_integer_, NA_integer_, 0L))

    out <- matchBarcodes(c("AAAAAA", "AAAGAA", "AAGAAA"), choices, reverse=TRUE, substitutions=2L)
    expect_identical(out$index, c(4L, 4L, 4L))
    expect_identical(out$mismatches, c(1L, 2L, 0L))
})
