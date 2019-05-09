# This tests the combineComboCounts function.
# library(testthat); library(gp.sa.screen); source("test-combine.R")

GENERATOR <- function(nsamples, npool, nreads) {
    output <- vector("list", nsamples)
    for (i in seq_along(output)) {
        X <- sample(npool, nreads, replace=TRUE)
        Y <- sample(npool, nreads, replace=TRUE)

        # Unique-ifying
        key <- paste0(X, ".", Y)
        tab <- table(key)
        re.X <- as.integer(sub("\\..*", "", names(tab)))
        re.Y <- as.integer(sub(".*\\.", "", names(tab)))

        output[[i]] <- DataFrame(combination=I(DataFrame(X=re.X, Y=re.Y)), count=as.integer(tab))
    }
    output
}

MANUAL <- function(...) {
    everything <- list(...)
    
    has.keys <- vector("list", length(everything))
    for (i in seq_along(has.keys)) {
        current <- everything[[i]]$combination
        has.keys[[i]] <- paste0(current[,1], ".", current[,2])
    }

    all.keys <- unique(unlist(has.keys))
    mat <- matrix(0L, length(all.keys), length(everything))
    for (i in seq_along(has.keys)) {
        m <- match(has.keys[[i]], all.keys)
        counts <- everything[[i]]$count
        mat[m,i] <- counts
    }

    out <- DataFrame(X=as.integer(sub("\\..*", "", all.keys)),
        Y=as.integer(sub(".*\\.", "", all.keys)))
    o <- order(out$X, out$Y)
    list(keys=out[o,], counts=mat[o,,drop=FALSE])
}

test_that("combineComboCounts works as expected", {
    for (N in c(1, 5, 10, 50, 100, 500)) {
        output <- GENERATOR(3, 10, N)
        out <- combineComboCounts(output[[1]], output[[2]], output[[3]])
        ref <- MANUAL(output[[1]], output[[2]], output[[3]])
        expect_equivalent(out$keys, ref$keys)
        expect_identical(out$counts, ref$counts)
    }
})

test_that("combineComboCounts works as expected", {
    N <- 500
    output <- GENERATOR(3, 10, N)
    out <- combineComboCounts(A=output[[1]], B=output[[2]], C=output[[3]])
    expect_identical(colnames(out$counts), c("A", "B", "C"))
})
