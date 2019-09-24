# Tests the combineBarcodeTests and barcodeSetTest functions.
# library(testthat); library(gp.sa.screen); source("setup.R"); source("test-consolidate.R")

##################################
##################################

stats <- DataFrame(PValue=runif(100), LogCPM=rnorm(100), LogFC=rnorm(100))
gene <- sample(c(letters, LETTERS), nrow(stats), replace=TRUE)

test_that("combineBarcodeTests works correctly", {
    per.gene <- combineBarcodeTests(stats, gene)
    expect_identical(rownames(per.gene), sort(unique(gene)))
    expect_true(all(c("PValue", "LogCPM", "LogFC", "FDR") %in% colnames(per.gene)))
    expect_true("Direction" %in% colnames(per.gene))
})

test_that("combineBarcodeTests handles NA's in the data.frame correctly", {
    old <- combineBarcodeTests(stats, gene)

    stats$PValue[1] <- NA
    per.gene <- combineBarcodeTests(stats, gene)
    expect_identical(rownames(per.gene), sort(unique(gene)))

    chosen <- gene[1]
    expect_identical(old[chosen,"NBarcodes"], per.gene[chosen,"NBarcodes"]+1L)

    stats$PValue[gene==chosen] <- NA
    per.gene <- combineBarcodeTests(stats, gene)
    expect_identical(0L, per.gene[chosen,"NBarcodes"])
})

test_that("combineBarcodeTests handles NA's in the genes correctly", {
    old <- combineBarcodeTests(stats, gene)

    discarded <- LETTERS[1:5]
    gene[gene %in% discarded] <- NA
    per.gene <- combineBarcodeTests(stats, gene)
    expect_identical(rownames(per.gene), sort(unique(setdiff(gene, discarded))))

    expect_identical(per.gene[,"PValue"], old[setdiff(rownames(old), discarded),"PValue"])
})

test_that("combineBarcodeTests handles multiple LogFCs correctly", {
    astats <- DataFrame(PValue=runif(100), LogCPM=rnorm(100), LogFC.1=rnorm(100), LogFC.2=rnorm(100))
    per.gene <- combineBarcodeTests(astats, gene)
    expect_identical(rownames(per.gene), sort(unique(gene)))
    expect_true(all(c("LogFC.1", "LogFC.2") %in% colnames(per.gene)))
    expect_false("Direction" %in% colnames(per.gene))
})

test_that("combineBarcodeTests works with the Holm-based approach correctly", {
    # Same as taking the minimum bonferroni p-value.
    gstats <- combineBarcodeTests(stats, gene, method="holm-min", min.sig.n=1, min.sig.prop=0)
    by.gene <- split(stats$PValue, gene)
    min.p <- vapply(by.gene, min, 0)
    expect_identical(unname(gstats$PValue), pmin(1, min.p * lengths(by.gene)))

    # Same as taking the maximum p-value.
    gstats <- combineBarcodeTests(stats, gene, method="holm-min", min.sig.n=1, min.sig.prop=1)
    max.p <- vapply(by.gene, FUN=function(p) max(p.adjust(p, "holm")), 0)
    expect_identical(gstats$PValue, max.p)

    gstats2 <- combineBarcodeTests(stats, gene, method="holm-min", min.sig.n=100, min.sig.prop=100)
    expect_identical(gstats, gstats2)

    # Same as taking the "middle" p-value.
    gstats <- combineBarcodeTests(stats, gene, method="holm-min", min.sig.n=1, min.sig.prop=0.5)
    mid.p <- vapply(by.gene, FUN=function(p) {
        p <- p.adjust(p, method="holm")
        if (length(p) %% 2) median(p) else sort(p)[length(p)/2]
    }, 0)
    expect_identical(gstats$PValue, mid.p)

    # Works correctly after shuffling the order.
    o <- sample(nrow(stats))
    gstats <- combineBarcodeTests(stats, gene, method="holm-min")
    gstats2 <- combineBarcodeTests(stats[o,], gene[o], method="holm-min")
    expect_identical(gstats, gstats2[rownames(gstats),])
})

test_that("combineBarcodeTests handles empty inputs correctly", {
    per.gene <- combineBarcodeTests(stats[0,], gene[0])
    expect_identical(nrow(per.gene), 0L)
})

##################################
##################################

N <- 100
mu <- 2^runif(N, 2, 10)
y <- matrix(rnbinom(N * 6, mu=mu, size=10), ncol=6)
g <- gl(2, 3)
design <- model.matrix(~g)
v <- limma::voom(y, design=design)

gene <- sample(LETTERS, N, replace=TRUE)
results <- DataFrame(LogFC=rnorm(N), LogCPM=runif(N))

test_that("barcodeSetTest works correctly", {
    per.gene <- barcodeSetTest(v, gene, design=design, contrast=2)
    expect_identical(rownames(per.gene), sort(unique(gene)))
    expect_true(all(c("PValue", "FDR") %in% colnames(per.gene)))
    expect_true("Direction" %in% colnames(per.gene))
})

test_that("barcodeSetTest works with results", {
    per.gene <- barcodeSetTest(v, gene, design=design, contrast=2, stats=results)
    expect_identical(rownames(per.gene), sort(unique(gene)))
    expect_true(all(c("PValue", "LogCPM", "LogFC", "FDR") %in% colnames(per.gene)))
    expect_true("Direction" %in% colnames(per.gene))

    expect_identical(per.gene$LogFC[1], median(results$LogFC[rownames(per.gene)[1]==gene]))
    expect_identical(per.gene$LogCPM[1], median(results$LogCPM[rownames(per.gene)[1]==gene]))
})

test_that("barcodeSetTest works with subsets", {
    chosen <- sample(N, 20)
    per.gene <- barcodeSetTest(v[chosen,], gene, design=design, contrast=2, stats=results, subset=chosen)
    expect_identical(rownames(per.gene), sort(unique(gene)))
    expect_true(all(c("PValue", "LogCPM", "LogFC", "FDR") %in% colnames(per.gene)))
    expect_true("Direction" %in% colnames(per.gene))

    ref <- barcodeSetTest(v[chosen,], gene[chosen], design=design, contrast=2, stats=results[chosen,])
    expect_identical(ref, per.gene[rownames(ref),])
    expect_true(all(is.na(per.gene[setdiff(rownames(per.gene), rownames(ref)),"PValue"])))
})
