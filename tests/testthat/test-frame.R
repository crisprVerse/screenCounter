# This tests the DiffScreenStatFrame class.
# library(testthat); library(gp.sa.screen); source("test-frame.R")

library(gp.sa.core)

test_that("constructors works as expected", {
    da.output <- DataFrame(LogFC=1:10, AveAb=1:10, FDR=0, PValue=0:9/10)
    Y <- DiffScreenStatFrame(da.output, 
        design=cbind(A=c(X=1, Y=-1), B=1),
        contrast=c(A=1, B=-1),
        method="voom", description="I did voom")

    expect_s4_class(Y, "DiffScreenStatFrame")
    expect_identical(trackinfo(Y)$method, "voom")
    expect_identical(trackinfo(Y)$description, "I did voom")
})

test_that(".trackCheck specialization works as expected for DiffScreenStatFrames", {
    tmp <- DataFrame(LogFC=1, AveAb=2, PValue=3, FDR=4)
    X <- as(tmp, "DiffScreenStatFrame")
    expect_error(.trackCheck(X), "origin")

    trackinfo(X) <- list(origin="SOME_RANDOM_ID", description="blah blah blah")
    expect_error(.trackCheck(X), "design")

    trackinfo(X)$design <- cbind(A=c(X=1, Y=-1), B=1)
    expect_error(.trackCheck(X), "contrast")

    trackinfo(X)$contrast <- c(A=1, B=-1)
    expect_error(.trackCheck(X), "method")

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
