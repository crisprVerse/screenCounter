# Tests the basic hashing capabilities for barcode searching.
# library(gp.sa.screen); library(testthat); source("test-hash.R")

BASES <- c("A", "C", "G", "T")

GENERATE_RANDOM_SEQ <- function(n) {
    paste(sample(BASES, n, replace=TRUE), collapse="")
}

MANUAL_HASH <- function(barcode, split=FALSE) {
    if (split) {
        B <- barcode
    } else { 
        B <- strsplit(barcode, "")[[1]]
    }

    collected <- numeric(0)
    while (length(B)) {
        word <- head(B, 16)
        as_num <- c(A=0, C=1, G=2, T=3)[word] * 4^(seq_along(word)-1)
        collected <- c(collected, sum(as_num))
        B <- tail(B, -16)
    }
    collected
}

test_that("basic hashing works correctly", {
    BASIC_COMPARE <- function(barcode) {
        out <- gp.sa.screen:::basic_hash(barcode)
        ref <- MANUAL_HASH(barcode)
        expect_identical(out, ref)
    }

    for (n in 1:40) {
        barcode <- GENERATE_RANDOM_SEQ(n)
        BASIC_COMPARE(barcode)
    }
})

test_that("shifted hashing works correctly", {
    SHIFT_COMPARE <- function(barcode, ensuing) {
        out <- gp.sa.screen:::shift_hash(barcode, ensuing)
        
        ensuing <- strsplit(ensuing, "")[[1]]
        collected <- vector("list", length(ensuing))
        for (i in seq_along(ensuing)) {
            barcode <- paste0(substring(barcode, 2, nchar(barcode)), ensuing[[i]])
            collected[[i]] <- MANUAL_HASH(barcode)
        }

        expect_identical(collected, out)
    }

    for (n1 in 1:40) {
        for (n2 in c(1, 5, 20, 40)) {
            barcode <- GENERATE_RANDOM_SEQ(n1)
            ensuing <- GENERATE_RANDOM_SEQ(n2)
            SHIFT_COMPARE(barcode, ensuing)
        }
    }
})

test_that("substituted hashing works correctly", {
    SUBSTITUTE_COMPARE <- function(barcode) {
        B <- strsplit(barcode, "")[[1]]

        collected <- list()
        for (pos in seq_along(B)) {
            for (alt in setdiff(BASES, B[[pos]])) {
                B2 <- B
                B2[[pos]] <- alt
                collected <- c(collected, list(MANUAL_HASH(B2, split=TRUE)))
            }
        }

        out <- gp.sa.screen:::substitute_hash(barcode)
        expect_identical(collected, out)
    }

    for (n in 1:40) {
        barcode <- GENERATE_RANDOM_SEQ(n)
        SUBSTITUTE_COMPARE(barcode)
    }
})

test_that("deleted hashing works correctly", {
    DELETE_COMPARE <- function(barcode) {
        B <- strsplit(barcode, "")[[1]]

        collected <- list()
        for (pos in rev(seq_along(B))) {
            B2 <- B[-pos]
            collected <- c(collected, list(MANUAL_HASH(B2, split=TRUE)))
        }

        out <- gp.sa.screen:::delete_hash(barcode)
        expect_identical(collected, out)
    }

    for (n in 1:40) {
        barcode <- GENERATE_RANDOM_SEQ(n)
        DELETE_COMPARE(barcode)
    }
})
