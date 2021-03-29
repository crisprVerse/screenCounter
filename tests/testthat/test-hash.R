# Tests the basic hashing capabilities for barcode searching.
# library(screenCounter); library(testthat); source("setup.R"); source("test-hash.R")

test_that("basic hashing works correctly", {
    for (n in 1:70) {
        barcode <- GENERATE_RANDOM_SEQ(n)
        out <- screenCounter:::basic_hash(barcode)
        ref <- MANUAL_HASH(barcode)
        expect_identical(out, ref)
    }
})

test_that("shifted hashing works correctly", {
    for (n1 in 1:70) {
        for (n2 in c(1, 5, 20, 40)) {
            barcode <- GENERATE_RANDOM_SEQ(n1)
            ensuing <- GENERATE_RANDOM_SEQ(n2)

            out <- screenCounter:::shift_hash(barcode, ensuing)
            
            ensuing <- strsplit(ensuing, "")[[1]]
            collected <- vector("list", length(ensuing))
            for (i in seq_along(ensuing)) {
                barcode <- paste0(substring(barcode, 2, nchar(barcode)), ensuing[[i]])
                collected[[i]] <- MANUAL_HASH(barcode)
            }

            expect_identical(collected, out)
        }
    }
})
