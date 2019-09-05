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
        commit="never", save.all=FALSE)

    expect_identical(names(out$barcode), NAME)
    expect_true(all(c("PValue", "FDR", "LogCPM", "LogFC") %in% colnames(out$barcode[[NAME]])))
    expect_identical(rownames(out$barcode[[NAME]]), NULL) # no names in 'se'.

    # All proper settings.
    out2 <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run",
        reference.field="time", reference.level=0,
        norm.type.field="class", norm.type.level="NEG",
        gene.field="gene",
        commit="never", save.all=FALSE
    )
    expect_identical(names(out2$barcode), NAME)
#    expect_identical(colnames(out2$barcode[[1]])[-1], colnames(out$barcode[[1]])) # skipping auto-added gene field.
    expect_identical(colnames(out2$barcode[[1]]), colnames(out$barcode[[1]])) 
    expect_identical(names(out2$gene), NAME)

    expect_true("PValue" %in% colnames(out2$gene[[NAME]]))
    expect_true("FDR" %in% colnames(out2$gene[[NAME]]))
    expect_true("LogCPM" %in% colnames(out2$gene[[NAME]]))
    expect_true("LogFC" %in% colnames(out2$gene[[NAME]]))
    expect_identical(rownames(out2$gene[[NAME]]), sort(unique(rowData(se)$gene)))

    # All default settings
    out3 <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run",
        reference.field=NULL, norm.type.field=NULL, gene.field="gene",
        commit="never", save.all=FALSE)
    expect_identical(names(out3), names(out2))
    expect_identical(colnames(out3$gene), colnames(out2$gene))
    expect_identical(colnames(out3$barcode), colnames(out2$barcode))
})

test_that("runVoomScreen works correctly with expansion of per-gene results", {
    se2 <- se
    discarded <- LETTERS[1:3]
    assay(se2)[rowData(se2)$gene %in% discarded,] <- 0

    out <- runVoomScreen(se2, covariates="time", comparisons=list("time"), block="run",
        reference.field="time", reference.level=0,
        norm.type.field="class", norm.type.level="NEG",
        gene.field="gene",
        commit="never", save.all=FALSE
    )

    expect_identical(names(out$barcode), NAME)
    expect_identical(names(out$gene), NAME)

    expect_identical(rownames(out$gene[[1]]), sort(unique(rowData(se)$gene)))
    expect_true(all(is.na(out$gene[[NAME]][discarded, "PValue"])))
})

test_that("runVoomScreen works correctly with other options", {
    ref <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run",
        reference.field="time", reference.level=0,
        norm.type.field="class", norm.type.level="NEG",
        gene.field="gene",
        commit="never", save.all=FALSE
    )

    # Checking that the annotation is properly stored.
    rowData(se)$SEVERANCE <- rowData(se)$gene
    alt <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run", lfc=1,
        reference.field="time", reference.level=0,
        norm.type.field="class", norm.type.level="NEG",
        gene.field="gene", annotation="SEVERANCE",
        commit="never", save.all=FALSE
    )
    expect_identical(alt[[1]][[1]]$SEVERANCE, rowData(se)$gene)
    expect_identical(alt[[2]][[1]]$SEVERANCE, rownames(alt[[2]][[1]]))

    # Checking that it responds to the lfc threshold.
    alt <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run", lfc=1,
        reference.field="time", reference.level=0,
        norm.type.field="class", norm.type.level="NEG",
        gene.field="gene",
        commit="never", save.all=FALSE
    )
    expect_true(isTRUE(all.equal(alt[[1]][[1]]$LogFC, ref[[1]][[1]]$LogFC)))
    expect_false(isTRUE(all.equal(alt[[1]][[1]]$PValue, ref[[1]][[1]]$PValue)))

    # Checking that it responds to contrasts.fun.
    alt <- runVoomScreen(se, covariates="time", block="run", 
        reference.field="time", reference.level=0,
        norm.type.field="class", norm.type.level="NEG",
        gene.field="gene", contrasts.fun=list(TIME=c(time=1)),
        commit="never", save.all=FALSE
    )
    expect_identical(names(alt[[1]]), "TIME")
    expect_identical(as.data.frame(alt[[1]][[1]]), as.data.frame(ref[[1]][[1]]))
    expect_identical(as.data.frame(alt[[2]][[1]]), as.data.frame(ref[[2]][[1]]))
})

test_that("runVoomScreen saves content correctly", {
    library(gp.sa.core)
    proj <- tempfile()
    newProject(proj, title="Aaron's project",
        description="This is one of Aaron's projects")

    report <- file.path(proj, "report.Rmd")
    trackinfo(se)$origin <- list(list(id="asdad"))

    out <- runVoomScreen(se, covariates="time", comparisons=list("time"), block="run",
        reference.field=NA, norm.type.field=NA, gene.field=NA,
        fname=report, save.all=TRUE)
    expect_true(file.exists(report))

    res.dir <- file.path(proj, "report-results")
    expect_true(file.exists(file.path(res.dir, "barcode.results-1")))
    expect_true(file.exists(getResultManifest(dir=res.dir)))

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
    expect_true(file.exists(getResultManifest(dir=res.dir)))
})
