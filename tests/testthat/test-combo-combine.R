# This tests the combineComboCounts function.
# library(testthat); library(screenCounter); source("test-combo-combine.R")

GENERATOR <- function(nsamples, npool, nreads) {
    output <- vector("list", nsamples)
    for (i in seq_along(output)) {
        X <- sample(npool, nreads, replace=TRUE)
        Y <- sample(npool, nreads, replace=TRUE)

        # Unique-ifying
        key <- paste0(X, ".", Y)
        tab <- table(key)

        re.X <- sub("\\..*", "", names(tab))
        re.Y <- sub(".*\\.", "", names(tab))
        if (is.integer(npool)) {
            re.X <- as.integer(re.X)
            re.Y <- as.integer(re.Y)
        }

        output[[i]] <- DataFrame(combinations=I(DataFrame(X=re.X, Y=re.Y)), counts=as.integer(tab))
    }
    output
}

MANUAL <- function(...) {
    everything <- list(...)
    
    has.keys <- vector("list", length(everything))
    for (i in seq_along(has.keys)) {
        current <- everything[[i]]$combinations
        has.keys[[i]] <- paste0(current[,1], ".", current[,2])
    }

    all.keys <- unique(unlist(has.keys))
    mat <- matrix(0L, length(all.keys), length(everything))
    for (i in seq_along(has.keys)) {
        m <- S4Vectors::match(has.keys[[i]], all.keys)
        counts <- everything[[i]]$counts
        mat[m,i] <- counts
    }

    out <- DataFrame(X=sub("\\..*", "", all.keys), Y=sub(".*\\.", "", all.keys))
    if (is.integer(everything[[1]]$combinations[,1])) {
        out[,1] <- as.integer(out[,1])
        out[,2] <- as.integer(out[,2])
    }

    o <- order(out$X, out$Y)
    list(combinations=out[o,], counts=mat[o,,drop=FALSE])
}

test_that("combineComboCounts works as expected with integers", {
    for (N in c(1, 5, 10, 50, 100, 500)) {
        output <- GENERATOR(3, 10, N)
        out <- combineComboCounts(output[[1]], output[[2]], output[[3]])
        ref <- MANUAL(output[[1]], output[[2]], output[[3]])
        expect_equivalent(out$combinations, ref$combinations)
        expect_identical(out$counts, ref$counts)
    }
})

test_that("combineComboCounts works as expected with strings", {
    for (N in c(1, 5, 10, 50, 100, 500)) {
        output <- GENERATOR(3, head(LETTERS, 10), N)
        out <- combineComboCounts(output[[1]], output[[2]], output[[3]])
        ref <- MANUAL(output[[1]], output[[2]], output[[3]])
        expect_equivalent(out$combinations, ref$combinations)
        expect_identical(out$counts, ref$counts)
    }
})

test_that("combineComboCounts respects file names", {
    N <- 500
    output <- GENERATOR(3, 10, N)
    out <- combineComboCounts(A=output[[1]], B=output[[2]], C=output[[3]])
    expect_identical(colnames(out$counts), c("A", "B", "C"))
})

test_that("combineComboCounts fails for non-unique entries", {
    N <- 500
    output <- GENERATOR(3, 10, N)
    expect_error(combineComboCounts(A=output[[1]][c(1,1:nrow(output[[1]])),], B=output[[2]], C=output[[3]]), "unique")
})
