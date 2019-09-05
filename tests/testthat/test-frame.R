# This tests the DAScreenStatFrame class.
# library(testthat); library(gp.sa.screen); source("test-frame.R")

# Mocking up an input into runVoom()
library(SummarizedExperiment)
se.input <- SummarizedExperiment()

library(gp.sa.core)
trackinfo(se.input)$origin <- list(list(id="SOME_ID"))

test_that("constructors works as expected", {
    da.output <- DataFrame(LogFC=1:10, LogCPM=1:10, FDR=0, PValue=0:9/10)
    Y <- DAScreenStatFrame(da.output, se.input,
        contrast=c(A=1, B=-1),
        method="voom", description="I did voom")

    expect_s4_class(Y, "DAScreenStatFrame")
    expect_identical(trackinfo(Y)$method, "voom")
    expect_identical(trackinfo(Y)$description, "I did voom")
    expect_identical(trackinfo(Y)$origin[[1]]$id, "SOME_ID")
})

test_that(".trackCheck specialization works as expected for DAScreenStatFrames", {
    tmp <- DataFrame(LogFC=1, LogCPM=2, PValue=3, FDR=4)
    X <- as(tmp, "DAScreenStatFrame")
    expect_error(.trackCheck(X), "origin")

    trackinfo(X) <- list(origin=list(list(id="SOME_RANDOM_ID")),
        description="blah blah blah")
    trackinfo(X)$contrast <- c(A=1, B=-1)
    expect_error(.trackCheck(X), "method")

    # Okay, adding a method.
    trackinfo(X)$method <- "random"
    expect_error(.trackCheck(X), "should be")

    # Okay, adding the right method.
    trackinfo(X)$method <- "voom"
    expect_error(.trackCheck(X), "feature")

    # Okay, adding a feature.
    trackinfo(X)$feature <- "blah"
    expect_error(.trackCheck(X), "should be either")

    trackinfo(X)$feature <- "gene"

    info <- .trackCheck(X)
    expect_identical(info$method, "voom")
    expect_match(info$type, "abundance")
})
