# Tests the dictionary building capabilities for barcode searching.
# library(screenCounter); library(testthat); source("setup.R"); source("test-build.R")

MANUAL_DICT <- function(barcodes, sub=FALSE, del=FALSE) {
    # Doesn't do anything special about clashes,
    # so make sure 'barcodes' doesn't have any clashes!
    if (!del) {
        idx <- seq_along(barcodes)
        priority <- integer(length(barcodes))
        keys <- lapply(barcodes, MANUAL_HASH)

        if (sub) {
            skeys <- lapply(barcodes, screenCounter:::substitute_hash)
            sidx <- rep(seq_along(skeys), lengths(skeys))

            keys <- c(keys, unlist(skeys, recursive=FALSE))
            idx <- c(idx, sidx)
            priority <- c(priority, rep(1L, length(sidx)))
        }
    } else {
        keys <- lapply(barcodes, screenCounter:::delete_hash)
        idx <- rep(seq_along(keys), lengths(keys))
        keys <- unlist(keys, recursive=FALSE)
        priority <- rep(1L, length(idx))
    }

    list(keys, list(idx, priority))
}

COMPARE_DICT <- function(left, right) {
    lvect <- do.call(mapply, c(left[[1]], list(FUN=c, SIMPLIFY=FALSE)))
    ol <- do.call(order, lvect)

    rvect <- do.call(mapply, c(right[[1]], list(FUN=c, SIMPLIFY=FALSE)))
    or <- do.call(order, rvect)
    
    # Comparing keys.
    expect_identical(
        lapply(lvect, "[", ol),
        lapply(rvect, "[", or)
    )

    expect_identical(left[[2]][[1]][ol], right[[2]][[1]][or]) # index
    expect_identical(left[[2]][[2]][ol], right[[2]][[2]][or]) # priority
}

test_that("dictionary building works correctly in the basic case", {
    for (n in c(10, 20, 40)) { 
        barcodes <- strrep(BASES, n) # avoid error from overlap.
        ref <- MANUAL_DICT(barcodes)
        out <- screenCounter:::build_dict(barcodes, FALSE, FALSE)
        expect_identical(as.integer(n), out[[2]])
        COMPARE_DICT(ref, out[[1]])
    }

    # Clashes cause errors.
    expect_error(screenCounter:::build_dict(c("AAAA", "AAAA"), FALSE, FALSE), "duplicated")

    # Variable length causes errors.
    expect_error(screenCounter:::build_dict(c("AAA", "AAAA"), FALSE, FALSE), "variable length")
})

test_that("dictionary building works correctly with substitutions", {
    for (n in c(10, 20, 40)) { 
        barcodes <- strrep(BASES, n)
        ref <- MANUAL_DICT(barcodes, sub=TRUE)
        out <- screenCounter:::build_dict(barcodes, TRUE, FALSE)
        expect_identical(as.integer(n), out[[2]])
        COMPARE_DICT(ref, out[[1]])
    }

    # Clashes between mismatch and perfect sequences from different barcodes are handled properly.
    truseq <- c("ACGTA", "ACGTG")
    out <- screenCounter:::build_dict(truseq, TRUE, FALSE)
    expect_identical(length(out[[1]][[1]]), sum(1L + 3L * nchar(truseq)) - 4L)

    for (i in seq_along(truseq)) {
        current <- MANUAL_HASH(truseq[i])
        idx <- which(vapply(out[[1]][[1]], FUN=function(x) identical(x, current), FUN.VALUE=TRUE))

        expect_identical(out[[1]][[2]][[1]][idx], i)
        expect_identical(out[[1]][[2]][[2]][idx], 0L)
    }

    # Clashes between two mismatch sequences are handled properly.
    truseq <- c("AAAAAACA", "AAAAAAAC")
    out <- screenCounter:::build_dict(truseq, TRUE, FALSE)

    for (conflict in c("AAAAAAAA", "AAAAAACC")) {
        conflictor <- MANUAL_HASH(conflict)
        chosen <- which(vapply(out[[1]][[1]], FUN=function(x) identical(x, conflictor), FUN.VALUE=TRUE))

        expect_identical(out[[1]][[2]][[1]][chosen], 0L) # Stored as -1 in C++, then +1 when moving back into R.
        expect_identical(out[[1]][[2]][[2]][chosen], 1L)
    }
})

test_that("dictionary building works correctly with deletions", {
    for (n in c(2, 5, 10)) { 
        # Carefully chosen to avoid consecutive runs that 
        # would result in two deletions clashing with each other.
        barcodes <- strrep(c("ACGT", "GCAT", "TACA", "GCTC"), n)
        ref <- MANUAL_DICT(barcodes, del=TRUE)
        out <- screenCounter:::build_dict(barcodes, FALSE, TRUE)
        expect_identical(unique(nchar(barcodes))-1L, out[[2]])
        COMPARE_DICT(ref, out[[1]])
    }

    # Clashes between two mismatch sequences are handled properly.
    truseq <- c("AAAAAACA", "AAAAAAAC")
    out <- screenCounter:::build_dict(truseq, FALSE, TRUE)

    for (conflict in c("AAAAAAA", "AAAAAAC")) {
        conflictor <- MANUAL_HASH(conflict)
        chosen <- which(vapply(out[[1]][[1]], FUN=function(x) identical(x, conflictor), FUN.VALUE=TRUE))
        expect_identical(out[[1]][[2]][[1]][chosen], 0L) # Stored as -1 in C++, then +1 when moving back into R.
        expect_identical(out[[1]][[2]][[2]][chosen], 1L)
    }

    # Clashes between two mismatch sequences of the same barcode do not create '-1's.
    truseq <- c("AAAAAAAA", "TTTTTTTT")
    out <- screenCounter:::build_dict(truseq, FALSE, TRUE)

    expect_identical(out[[1]][[2]][[2]], rep(1L, length(truseq)))
    out[[1]][[2]][[2]][] <- 0L
    COMPARE_DICT(out[[1]], MANUAL_DICT(c("AAAAAAA", "TTTTTTT")))

    # Variable length causes errors.
    expect_error(screenCounter:::build_dict(c("AAA", "AAAA"), FALSE, TRUE), "variable length")
})
