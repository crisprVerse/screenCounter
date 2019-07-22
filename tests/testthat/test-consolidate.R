# Tests the barcode2genes function.
# library(testthat); library(gp.sa.screen); source("setup.R"); source("test-consolidate.R")

stats <- DataFrame(PValue=runif(100), LogCPM=rnorm(100), LogFC=rnorm(100))
gene <- sample(c(letters, LETTERS), nrow(stats), replace=TRUE)

test_that("barcode2genes works correctly", {
    per.gene <- barcodes2genes(stats, gene)                  
    expect_identical(rownames(per.gene), sort(unique(gene)))
    expect_true(all(c("PValue", "LogCPM", "LogFC", "FDR") %in% colnames(per.gene)))
    expect_true("Direction" %in% colnames(per.gene))
})

test_that("barcode2genes handles NA's in the data.frame correctly", {
    old <- barcodes2genes(stats, gene)

    stats$PValue[1] <- NA
    per.gene <- barcodes2genes(stats, gene)
    expect_identical(rownames(per.gene), sort(unique(gene)))

    chosen <- gene[1]
    expect_identical(old[chosen,"NBarcodes"], per.gene[chosen,"NBarcodes"]+1L)
})

test_that("barcode2genes handles NA's in the genes correctly", {
    old <- barcodes2genes(stats, gene)

    discarded <- LETTERS[1:5]
    gene[gene %in% discarded] <- NA
    per.gene <- barcodes2genes(stats, gene)
    expect_identical(rownames(per.gene), sort(unique(setdiff(gene, discarded))))

    expect_identical(per.gene[,"PValue"], old[setdiff(rownames(old), discarded),"PValue"])
})

test_that("barcode2genes handles multiple LogFCs correctly", {
    astats <- DataFrame(PValue=runif(100), LogCPM=rnorm(100), LogFC.1=rnorm(100), LogFC.2=rnorm(100))
    per.gene <- barcodes2genes(astats, gene)
    expect_identical(rownames(per.gene), sort(unique(gene)))
    expect_true(all(c("LogFC.1", "LogFC.2") %in% colnames(per.gene)))
    expect_false("Direction" %in% colnames(per.gene))
})

test_that("barcode2genes handles empty inputs correctly", {
    per.gene <- barcodes2genes(stats[0,], gene[0])
    expect_identical(nrow(per.gene), 0L)
})

