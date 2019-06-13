# This runs through all of the voom-specific functionality.
# library(testthat); library(gp.sa.screen); source("test-voom.R")

set.seed(1000)

N <- 10000
mu <- 2^runif(N, 3, 5)
group <- rep(0:2, each=3)

library(SummarizedExperiment)
se <- SummarizedExperiment(
    list(counts=matrix(rnbinom(N*length(group), mu=mu, size=10), nrow=N)), 
    colData=DataFrame(time=group, run=factor(rep(1:3, 3))),
    rowData=DataFrame(gene=sample(c(NA_character_, letters, LETTERS), N, replace=TRUE))
)
rowData(se)$class <- ifelse(is.na(rowData(se)$gene), "NEG", ".")

TB_NAME <- "Effect of increasing `time` (de) (barcode)"
TG_NAME <- "Effect of increasing `time` (de) (gene)"

test_that("runVoomScreen works correctly in basic scenarios", {
    # All missing settings.
    out <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run",
        reference.field=NA, norm.type.field=NA, gene.field=NA)
    expect_identical(names(out), TB_NAME)

    expect_true("PValue" %in% colnames(out[[TB_NAME]]))
    expect_true("FDR" %in% colnames(out[[TB_NAME]]))
    expect_true("AverageAbundance" %in% colnames(out[[TB_NAME]]))
    expect_true("LogFC" %in% colnames(out[[TB_NAME]]))
    expect_identical(rownames(out[[TB_NAME]]), NULL) # no names in 'se'.

    # All proper settings.
    out <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run",
        reference.field="time", reference.level=0,
        norm.type.field="class", norm.type.level="NEG",
        gene.field="gene"
    )
    expect_identical(names(out), c(TB_NAME, TG_NAME))

    expect_true("PValue" %in% colnames(out[[TG_NAME]]))
    expect_true("FDR" %in% colnames(out[[TG_NAME]]))
    expect_true("BestAverageAbundance" %in% colnames(out[[TG_NAME]]))
    expect_true("BestLogFC" %in% colnames(out[[TG_NAME]]))
    expect_identical(rownames(out[[TG_NAME]]), sort(unique(rowData(se)$gene)))

    # All default settings
    out2 <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run",
        reference.field=NULL, norm.type.field=NULL, gene.field="gene")
    expect_identical(names(out2), c(TB_NAME, TG_NAME))

    expect_identical(nrow(out2[[TB_NAME]]), nrow(out[[TB_NAME]]))
    expect_identical(rownames(out2[[TG_NAME]]), rownames(out[[TG_NAME]]))

    expect_false(isTRUE(all.equal(out2[[TB_NAME]], out[[TB_NAME]])))
    expect_false(isTRUE(all.equal(out2[[TG_NAME]], out[[TG_NAME]])))
})

test_that("runVoomScreen works correctly with expansion of per-gene results", {
    se2 <- se
    discarded <- LETTERS[1:3]
    assay(se2)[rowData(se2)$gene %in% discarded,] <- 0

    out <- runVoomScreen(se2, covariates="time", comparisons=list("time"), block="run",
        reference.field="time", reference.level=0,
        norm.type.field="class", norm.type.level="NEG",
        gene.field="gene"
    )
    expect_identical(names(out), c(TB_NAME, TG_NAME))

    expect_identical(rownames(out[[TG_NAME]]), sort(unique(rowData(se)$gene)))
    expect_true(all(is.na(out[[TG_NAME]][discarded, "PValue"])))
})


