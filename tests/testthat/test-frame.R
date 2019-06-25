# This tests the DBAStatFrame and DGAStatFrame class.
# library(testthat); library(gp.sa.screen); source("test-frame.R")

# Mocking up an input into runVoom()
library(SummarizedExperiment)
se.input <- SummarizedExperiment()

library(gp.sa.core)
trackinfo(se.input)$origin <- list(id="SOME_ID")

# Mocking up the aftermath of a DA analysis:
da.output <- DataFrame(N=1:10, PValue=0:9/10)

test_that("constructors works as expected", {
    Y <- DBAStatFrame(da.output, se.input,
        contrast=c(A=1, B=-1),
        method="voom", description="I did voom")

    expect_s4_class(Y, "DBAStatFrame")
    expect_identical(trackinfo(Y)$method, "voom")
    expect_identical(trackinfo(Y)$description, "I did voom")
    expect_identical(trackinfo(Y)$origin$id, "SOME_ID")

    Y <- DGAStatFrame(da.output, se.input,
        contrast=c(A=1, B=-1),
        method="voom", description="I did voom")

    expect_s4_class(Y, "DGAStatFrame")
    expect_identical(trackinfo(Y)$method, "voom")
    expect_identical(trackinfo(Y)$description, "I did voom")
    expect_identical(trackinfo(Y)$origin$id, "SOME_ID")
})

test_that("trackcheck specialization works as expected for DGAStatFrames", {
    tmp <- DataFrame(a=1, b=2, c="3")
    X <- as(tmp, "DGAStatFrame")
    expect_error(trackcheck(X), "origin")

    trackinfo(X) <- list(origin=list(id="SOME_RANDOM_ID"),
        description="blah blah blah")
    trackinfo(X)$contrast <- c(A=1, B=-1)

    # Okay, adding a method.
    expect_error(trackcheck(X), "that was used")
    trackinfo(X)$method <- "random"
    expect_error(trackcheck(X), "should be 'voom'")

    # Fine!
    trackinfo(X)$method <- "voom"
    info <- trackcheck(X)
    expect_identical(info$method, "voom")
    expect_match(info$type, "gene abundance")
})

test_that("trackcheck specialization works as expected for DBAStatFrames", {
    tmp <- DataFrame(a=1, b=2, c="3")
    X <- as(tmp, "DBAStatFrame")
    expect_error(trackcheck(X), "origin")

    trackinfo(X) <- list(origin=list(id="SOME_RANDOM_ID"),
        description="blah blah blah")
    trackinfo(X)$contrast <- c(A=1, B=-1)

    # Okay, adding a method.
    trackinfo(X)$method <- "random"
    expect_error(trackcheck(X), "should be 'voom'")

    # Fine!
    trackinfo(X)$method <- "voom"
    info <- trackcheck(X)
    expect_identical(info$method, "voom")
    expect_match(info$type, "barcode abundance")
})
