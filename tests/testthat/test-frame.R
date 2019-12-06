# This tests the ScreenBarcodeStatFrame class.
# library(testthat); library(gp.sa.screen); source("test-frame.R")

library(gp.sa.core)

test_that("constructors works as expected", {
    da.output <- DataFrame(LogFC=1:10, AveAb=1:10, FDR=0, PValue=0:9/10)
    Y <- ScreenBarcodeStatFrame(da.output, 
        design=cbind(A=c(X=1, Y=-1), B=1),
        contrast=c(A=1, B=-1),
        method="voom", description="I did voom")

    expect_s4_class(Y, "ScreenBarcodeStatFrame")
    expect_identical(trackinfo(Y)$method, "voom")
    expect_identical(trackinfo(Y)$description, "I did voom")
})

test_that(".trackCheck specialization works as expected for ScreenBarcodeStatFrames", {
    tmp <- DataFrame(LogFC=1, AveAb=2, PValue=3, FDR=4)
    X <- as(tmp, "ScreenBarcodeStatFrame")
    expect_error(.trackCheck(X), "origin")

    trackinfo(X)$origin <- "SOME_RANDOM_ID"
    trackinfo(X)$description <- "blah blah blah"
    expect_error(.trackCheck(X), "design")

    trackinfo(X)$design <- cbind(A=c(X=1, Y=-1), B=1)
    expect_error(.trackCheck(X), "contrast")

    trackinfo(X)$contrast <- c(A=1, B=-1)
    expect_error(.trackCheck(X), "method")

    # Okay, adding the right method.
    trackinfo(X)$method <- "voom"

    info <- .trackCheck(X)
    expect_identical(info$method, "voom")
    expect_match(info$type, "abundance")
})

test_that(".trackCheck specialization works as expected for ScreenFeatureStatFrames", {
    tmp <- DataFrame(LogFC=1, AveAb=2, PValue=3, FDR=4)
    X <- as(tmp, "ScreenFeatureStatFrame")
    trackinfo(X)$origin <- "SOME_RANDOM_ID"
    trackinfo(X)$description <- "blah blah blah"
    trackinfo(X)$design <- cbind(A=c(X=1, Y=-1), B=1)
    trackinfo(X)$contrast <- c(A=1, B=-1)
    trackinfo(X)$method <- "voom"

    expect_error(.trackCheck(X), "NBarcodes")
    X$NBarcodes <- 1L
    expect_error(.trackCheck(X), NA)
})

test_that(".trackCheck checks the row names correctly", {
    old_d <- .getDimnames()

    tmp <- DataFrame(LogFC=1, AveAb=2, PValue=3, FDR=4)
    Y <- ScreenBarcodeStatFrame(tmp,
        design=cbind(A=c(X=1, Y=-1), B=1),
        contrast=c(A=1, B=-1), feature="barcode",
        origin="YAY", method="voom", description="I did voom")
    rownames(Y) <- "X"

    expect_warning(.trackCheck(Y), "non-permissible")

    .setDimnames(old_d)
})
