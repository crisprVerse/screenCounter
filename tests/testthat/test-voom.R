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

NAME <- "Effect of increasing `time`"

test_that("runVoomScreen works correctly in basic scenarios", {
    # All missing settings.
    out <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run",
        reference.field=NA, norm.type.field=NA, gene.field=NA, 
        commit="never")

    expect_identical(names(out$barcode), NAME)
    expect_true(all(c("PValue", "FDR", "AveAb", "LogFC") %in% colnames(out$barcode[[NAME]])))
    expect_identical(rownames(out$barcode[[NAME]]), NULL) # no names in 'se'.

    # All proper settings.
    out2 <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run",
        reference.field="time", reference.level=0,
        norm.type.field="class", norm.type.level="NEG",
        gene.field="gene",
        commit="never"
    )
    expect_identical(names(out2$barcode), NAME)
#    expect_identical(colnames(out2$barcode[[1]])[-1], colnames(out$barcode[[1]])) # skipping auto-added gene field.
    expect_identical(colnames(out2$barcode[[1]]), colnames(out$barcode[[1]])) 
    expect_identical(names(out2$gene), NAME)

    expect_true("PValue" %in% colnames(out2$gene[[NAME]]))
    expect_true("FDR" %in% colnames(out2$gene[[NAME]]))
    expect_true("AveAb" %in% colnames(out2$gene[[NAME]]))
    expect_true("LogFC" %in% colnames(out2$gene[[NAME]]))
    expect_identical(rownames(out2$gene[[NAME]]), sort(unique(rowData(se)$gene)))

    # All default settings
    out3 <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run",
        reference.field=NULL, norm.type.field=NULL, gene.field="gene",
        commit="never")
    expect_identical(names(out3), names(out2))
    expect_identical(colnames(out3$gene), colnames(out2$gene))
    expect_identical(colnames(out3$barcode), colnames(out2$barcode))

    # Using other consolidation strategies.
    out3b <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run",
        reference.field=NULL, norm.type.field=NULL, gene.field="gene", method="fry",
        commit="never")
    expect_identical(names(out3b), names(out2))
    expect_identical(colnames(out3b$gene), colnames(out2$gene))
    expect_identical(colnames(out3b$barcode), colnames(out2$barcode))
    expect_true(identical(out3$barcode, out3b$barcode))
    expect_false(identical(out3$gene, out3b$gene))
})

test_that("runVoomScreen works with subsetting", {
    se$keep <- rep(c(TRUE, TRUE, FALSE), 3)
    ref <- runVoomScreen(se[,se$keep], covariates="time", comparisons=list("time"), block="run",
        reference.field=NULL, norm.type.field=NULL, gene.field="gene",
        commit="never")
    out <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run",
        reference.field=NULL, norm.type.field=NULL, gene.field="gene",
        commit="never", subset.factor="keep", subset.levels=TRUE)

    expect_match(trackinfo(out[[1]][[1]])$subset, "keep")
    expect_match(trackinfo(out[[2]][[1]])$subset, "keep")
    trackinfo(out[[1]][[1]])$subset <- NULL
    trackinfo(out[[2]][[1]])$subset <- NULL
    expect_identical(ref, out)
})

test_that("runVoomScreen works correctly with expansion of per-gene results", {
    se2 <- se
    discarded <- LETTERS[1:3]
    assay(se2)[rowData(se2)$gene %in% discarded,] <- 0

    out <- runVoomScreen(se2, covariates="time", comparisons=list("time"), block="run",
        reference.field="time", reference.level=0,
        norm.type.field="class", norm.type.level="NEG",
        gene.field="gene",
        commit="never"
    )

    expect_identical(names(out$barcode), NAME)
    expect_identical(names(out$gene), NAME)
    expect_identical(rownames(out$gene[[1]]), sort(unique(rowData(se)$gene)))
    expect_true(all(is.na(out$gene[[NAME]][discarded, "PValue"])))

    out2 <- runVoomScreen(se2, covariates="time", comparisons=list("time"), block="run",
        reference.field="time", reference.level=0,
        norm.type.field="class", norm.type.level="NEG",
        gene.field="gene", method="fry",
        commit="never"
    )

    expect_identical(names(out2$barcode), NAME)
    expect_identical(names(out2$gene), NAME)
    expect_identical(sort(rownames(out2$gene[[1]])), sort(unique(rowData(se)$gene)))
    expect_true(all(is.na(out2$gene[[NAME]][discarded, "PValue"])))
})

test_that("runVoomScreen works correctly with other options", {
    ref <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run",
        reference.field="time", reference.level=0,
        norm.type.field="class", norm.type.level="NEG",
        gene.field="gene",
        commit="never"
    )

    # Checking that the annotation is properly stored.
    rowData(se)$SEVERANCE <- rowData(se)$gene
    alt <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run", lfc=1,
        reference.field="time", reference.level=0,
        norm.type.field="class", norm.type.level="NEG",
        gene.field="gene", annotation="SEVERANCE",
        commit="never"
    )
    expect_identical(alt[[1]][[1]]$SEVERANCE, rowData(se)$gene)
    expect_identical(alt[[2]][[1]]$SEVERANCE, rownames(alt[[2]][[1]]))

    # Checking that it responds to the lfc threshold.
    alt <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run", lfc=1,
        reference.field="time", reference.level=0,
        norm.type.field="class", norm.type.level="NEG",
        gene.field="gene",
        commit="never"
    )
    expect_true(isTRUE(all.equal(alt[[1]][[1]]$LogFC, ref[[1]][[1]]$LogFC)))
    expect_false(isTRUE(all.equal(alt[[1]][[1]]$PValue, ref[[1]][[1]]$PValue)))

    # Checking that it responds to contrasts=.
    alt <- runVoomScreen(se, covariates="time", block="run", 
        reference.field="time", reference.level=0,
        norm.type.field="class", norm.type.level="NEG",
        gene.field="gene", contrasts=list(TIME=c(time=1)),
        commit="never"
    )
    expect_identical(names(alt[[1]]), "TIME")
    expect_identical(as.data.frame(alt[[1]][[1]]), as.data.frame(ref[[1]][[1]]))
    expect_identical(as.data.frame(alt[[2]][[1]]), as.data.frame(ref[[2]][[1]]))

    # Checking that it uses a different consolidation strategy.
    alt <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run",
        reference.field="time", reference.level=0,
        norm.type.field="class", norm.type.level="NEG",
        gene.field="gene", method="holm-mid",
        commit="never"
    )
    expect_true(identical(alt$barcode[[1]], ref$barcode[[1]]))
    expect_false(identical(alt$gene[[1]], ref$gene[[1]]))
})

test_that("runVoomScreen saves content correctly", {
    library(gp.sa.core)
    proj <- tempfile()
    newProject(proj, title="Aaron's project",
        description="This is one of Aaron's projects")

    report <- file.path(proj, "report.Rmd")
    trackinfo(se)$origin <- "asdad"

    out <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run",
        reference.field=NA, norm.type.field=NA, gene.field=NA,
        fname=report, save.all=TRUE)
    expect_true(file.exists(report))

    res.dir <- file.path(proj, "report-results")
    expect_true(file.exists(file.path(res.dir, "barcode.results-1")))
    expect_true(length(readResultManifest(dir=res.dir))>0)

    # Trying to save with gene-level information now.
    out <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run",
        reference.field="time", reference.level=0,
        norm.type.field="class", norm.type.level="NEG",
        gene.field="gene",
        fname=report, commit="never"
    )

    res.dir <- file.path(proj, "report-results")
    expect_true(file.exists(file.path(res.dir, "barcode.results-1")))
    expect_true(file.exists(file.path(res.dir, "gene.results-1")))
    expect_true(length(readResultManifest(dir=res.dir))>0)

    # Trying to save with fry.
    out <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run",
        reference.field="time", reference.level=0,
        norm.type.field="class", norm.type.level="NEG",
        gene.field="gene", method="fry",
        fname=report, commit="never"
    )

    res.dir <- file.path(proj, "report-results")
    expect_true(file.exists(file.path(res.dir, "barcode.results-1")))
    expect_true(file.exists(file.path(res.dir, "gene.results-1")))
    expect_true(length(readResultManifest(dir=res.dir))>0)
   
    # Dumps out normalized expression values.
    out <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run",
        reference.field=NA, norm.type.field=NA, gene.field=NA,
        fname=report, dump.norm="whee.csv", commit="never")
    expect_true(file.exists("whee.csv"))
})
