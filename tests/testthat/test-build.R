# Tests the dictionary building capabilities for barcode searching.
# library(screenCounter); library(testthat); source("setup.R"); source("test-build.R")

test_that("basic dictionary construction works as expected", {
    for (n1 in 5*1:15) {
        barcode <- vapply(1:20, function(i) GENERATE_RANDOM_SEQ(n1), "")
        barcode <- unique(barcode)
        dict <- screenCounter:::export_basic_dictionary(barcode)

        hashes <- dict[[1]]
        indices <- dict[[2]]
        expect_identical(sort(indices), seq_along(indices))

        for (x in seq_along(hashes)) {
            idx <- indices[[x]]
            expect_identical(hashes[[x]], MANUAL_HASH(barcode[idx]))
        }
    }

    expect_error(screenCounter:::export_basic_dictionary(rep("ACAT", 2)), "duplicated")
})

test_that("combined dictionary construction works as expected (basic)", {
    constant <- c("ACGT", "TCGA")
    barcode <- c("ACACA", "GTGTG", "CTCTC")
    dict <- screenCounter:::export_combined_dictionary(constant, list(barcode), 0, 0, 0, 0)

    expect_identical(length(dict), 1L)
    expect_identical(dict[[1]][[1]], 13L) # single length
    expect_identical(dict[[1]][[3]], integer(3)) # no mismatches

    idx <- dict[[1]][[4]]
    expect_identical(lengths(idx), rep(1L, 3))
    expect_identical(sort(unlist(idx)), 1:3)

    full.seq <- paste0(constant[1], barcode, constant[2])
    expect_identical(dict[[1]][[2]], lapply(full.seq[unlist(idx)], MANUAL_HASH))

    expect_error(screenCounter:::export_combined_dictionary(constant, list(c(barcode, barcode)), 0, 0, 0, 0), "duplicated")
    expect_error(screenCounter:::export_combined_dictionary(constant, list(c("A", barcode)), 0, 0, 0, 0), "variable length")
})

test_that("combined dictionary construction works as expected (multiple)", {
    constant <- c("ACGT", "ATTA", "TCGA")
    barcode1 <- c("ACACA", "GTGTG", "CTCTC")
    barcode2 <- c("TCTC", "TGTG")
    dict <- screenCounter:::export_combined_dictionary(constant, list(barcode1, barcode2), 0, 0, 0, 0)

    expect_identical(length(dict), 1L)
    dict <- dict[[1]]
    expect_identical(dict[[1]], 21L)

    n1 <- length(barcode1)
    n2 <- length(barcode2)
    expect_identical(length(dict[[3]]), n1*n2)
    expect_identical(dict[[3]], integer(n1*n2)) # no mismatches

    idx <- do.call(rbind, dict[[4]])
    expect_identical(nrow(idx), n1*n2)
    o <- order(idx[,2], idx[,1])
    expect_identical(idx[o,], unname(as.matrix(expand.grid(seq_len(n1), seq_len(n2)))))

    full.seq <- paste0(constant[1], barcode1[idx[,1]], constant[2], barcode2[idx[,2]], constant[3])
    expect_identical(dict[[2]], lapply(full.seq, MANUAL_HASH))

    expect_error(screenCounter:::export_combined_dictionary(constant, list(barcode1, c(barcode2, barcode2)), 0, 0, 0, 0), "duplicated")
    expect_error(screenCounter:::export_combined_dictionary(constant, list(barcode1), 0, 0, 0, 0), "number of constant regions")
})

library(S4Vectors)
COMPARE_HASHES <- function(x, y) {
    x1 <- DataFrame(do.call(rbind, x))
    y1 <- DataFrame(do.call(rbind, y))
    expect_identical(sort(x1), sort(y1))
}

COMPUTE_NSUBS <- function(full, N) {
    length(full) * vapply(seq_len(N), function(n) choose(nchar(full)[1], n) * 3L^n, 0)
}

test_that("combined dictionary construction works as expected (single substitution)", {
    # NOTE: avoid choosing barcodes that causes an overlap on substitution.
    constant <- c("ACGT", "TCGA")
    barcode <- c("ACACA", "GTGTG", "CACAC")

    dict <- screenCounter:::export_combined_dictionary(constant, list(barcode), 1, 0, 0, 1)
    expect_identical(length(dict), 1L)

    # Mutate the sequence and check.
    all.full <- paste0(constant[1], barcode, constant[2])
    for (i in seq_along(barcode)) {
        full <- all.full[i]
        len <- nchar(full)

        current <- character()
        for (j in seq_len(len)) {
            base <- substr(full, j, j)
            current <- c(current, paste0(substr(full, 1, j-1), setdiff(c("A", "C", "G", "T"), base), substr(full, j+1, len)))
        }

        keep <- unlist(dict[[1]][[4]])==i 
        expect_identical(MANUAL_HASH(full), dict[[1]][[2]][keep & dict[[1]][[3]]==0][[1]]) 

        observed <- dict[[1]][[2]][keep & dict[[1]][[3]]==1]
        ref <- lapply(current, MANUAL_HASH)
        COMPARE_HASHES(observed, ref)
    }

    expect_equal(length(dict[[1]][[2]]), sum(COMPUTE_NSUBS(all.full, 1L)) + length(all.full))
})

test_that("combined dictionary construction works as expected (single deletions)", {
    # NOTE: avoid choosing barcodes that causes an overlap on deletion.
    constant <- c("ACGT", "TCGA")
    barcode <- c("ACACA", "GTGTG", "CGCGC")

    dict <- screenCounter:::export_combined_dictionary(constant, list(barcode), 0, 0, 1, 1)
    expect_identical(length(dict), 2L)

    # Mutate the sequence and check.
    all.full <- paste0(constant[1], barcode, constant[2])
    for (i in seq_along(barcode)) {
        full <- all.full[i]
        len <- nchar(full)

        current <- character()
        for (j in seq_len(len)) {
            current <- c(current, paste0(substr(full, 1, j-1), substr(full, j+1, len)))
        }

        keep <- unlist(dict[[1]][[4]]) == i 
        expect_true(all(dict[[1]][[3]][keep] == 1))
        observed <- dict[[1]][[2]][keep]
        ref <- lapply(current, MANUAL_HASH)
        COMPARE_HASHES(observed, ref)
    }

    # The original sequences work correctly.
    expect_identical(lapply(all.full[unlist(dict[[2]][[4]])], MANUAL_HASH), dict[[2]][[2]])
})

test_that("combined dictionary construction works as expected (single insertion)", {
    # NOTE: avoid choosing barcodes that causes an overlap on insertion.
    constant <- c("ACGT", "TCGA")
    barcode <- c("ACACA", "GTGTG", "CGCGC")

    dict <- screenCounter:::export_combined_dictionary(constant, list(barcode), 0, 1, 0, 1)
    expect_identical(length(dict), 2L)

    # Mutate the sequence and check.
    all.full <- paste0(constant[1], barcode, constant[2])
    for (i in seq_along(barcode)) {
        full <- all.full[i]
        len <- nchar(full)

        current <- character()
        for (j in seq_len(len)[-1]) { # exclude the position at the start.
            current <- c(current, paste0(substr(full, 1, j-1), c("A", "C", "G", "T"), substr(full, j, len)))
        }

        keep <- unlist(dict[[1]][[4]]) == i 
        expect_true(all(dict[[1]][[3]][keep] == 1))
        observed <- dict[[1]][[2]][keep]
        ref <- lapply(unique(current), MANUAL_HASH)
        COMPARE_HASHES(observed, ref)
    }

    # The original sequences work correctly.
    expect_identical(lapply(all.full[unlist(dict[[2]][[4]])], MANUAL_HASH), dict[[2]][[2]])
})

test_that("combined dictionary construction works as expected (multiple edits)", {
    # NOTE: 'barcode' must be carefully chosen as two substitutions can cause an overlap.
    constant <- c("ACGT", "TCGA")
    barcode <- c("ACACA", "GTGTG", "CACAC")
    dict <- screenCounter:::export_combined_dictionary(constant, list(barcode), 2, 0, 0, 2)
    dict <- dict[[1]]

    full <- paste0(constant[1], barcode, constant[2])
    n <- COMPUTE_NSUBS(full, 2L)
    expect_equal(sum(dict[[3]] == 0L), length(barcode))
    expect_equal(sum(dict[[3]] == 1L), n[1])
    expect_equal(sum(dict[[3]] == 2L), n[2])

    # Total edits has some effect.
    ref <- screenCounter:::export_combined_dictionary(constant, list(barcode), 1, 0, 0, 1)
    dict <- screenCounter:::export_combined_dictionary(constant, list(barcode), 2, 0, 0, 1)
    COMPARE_HASHES(ref[[1]][[2]], dict[[1]][[2]])

    # Throwing them all together seems to do something. 
    dict <- screenCounter:::export_combined_dictionary(constant, list(barcode), 2, 1, 1, 2)
    expect_true(length(dict) > length(ref))
    offender <- which(vapply(dict, function(x) x[[1]], 0L)==ref[[1]][[1]])
    expect_true(length(dict[[offender]][[2]]) > length(ref[[1]][[2]]))
})

test_that("combined dictionary construction works as expected (overlaps)", {
    dict <- screenCounter:::export_combined_dictionary(c("A", "A"), list("A"), 0, 1, 1, 2)

    # Range of elements available.
    spread <- vapply(dict, function(x) x[[1]], 0L)
    expect_identical(sort(spread), 2:4)

    # Exact match is present.
    offender <- which(spread==3L)
    exactor <- dict[[offender]]
    expect_true(any(exactor[[3]]==0))

    # Mimic substitutions at every base except the first (which insertions can't touch).
    other.bases <- c("C", "T", "G")
    others <- c(paste0("A", other.bases, "A"), paste0("AA", other.bases))
    COMPARE_HASHES(exactor[[2]][exactor[[3]]!=0], lapply(others, MANUAL_HASH))
})
