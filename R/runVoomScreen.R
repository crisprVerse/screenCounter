#' Use \code{voom} to analyze a screen
#'
#' Perform a \code{voom} analysis on a count matrix from a sequencing screen to detect differentially abundant barcodes across samples.
#'
#' @param se A \linkS4class{SummarizedExperiment} object containing read counts for each barcode (row) and sample (column).
#' Alternatively, a database ID string or list, see \code{?"\link{gp.sa.core-data-inputs}"}.
#' @param ... Further arguments to pass to \code{\link{runVoom}}. 
#' @param reference.field String specifying the column of \code{colData(se)} containing the type of each sample (i.e., reference or not).
#' @param reference.level Character vector specifying the reference levels of the column named by \code{reference.field}.
#' @param norm.type.field String specifying the field of \code{rowData(se)} containing the gene type for each barcode.
#' @param norm.type.level Character vector specifying the gene types on which to perform normalization.
#' @param gene.field String specifying the field of \code{rowData(se)} that contains the gene identifier for each barcode.
#' @param method String specifying the consolidation method to convert per-barcode statistics into per-gene results.
#' @param save.all Logical scalar indicating whether the returned \linkS4class{DAScreenStatFrame}s should also be saved to file.
#' Defaults to \code{FALSE} if \code{se} is a SummarizedExperiment without provenance information, and \code{TRUE} otherwise.
#' @param dump.norm String specifying a path to an output file to save normalized abundances in a CSV file.
#' This file is intended only for diagnostic inspection and should \emph{not} be used as input into further GPSA pipeline.
#' @inheritParams gp.sa.diff::runVoom
#'
#' @return A \linkS4class{List} containing two Lists, \code{barcode} and \code{gene}.
#' Each list contains barcode- and gene-level result tables as \linkS4class{DAScreenStatFrame}s from all contrasts.
#' If \code{gene.field=NA}, only the \code{barcode} List is returned.
#' 
#' A Rmarkdown file is also created at \code{fname}, containing the steps required to reproduce the analysis.
#' This also provides a basis for further customization.
#'
#' @details
#' This function is largely a wrapper around \code{\link{runVoom}}, with three additional features:
#' \itemize{
#' \item Filtering based on reference samples, see below.
#' If \code{reference.field=NULL}, default edgeR filtering is used instead based \code{\link{filterByExpr}}.
#' If \code{reference.field=NA}, no filtering is performed.
#' \item Normalization based on non-targeting genes (NTGs) and/or non-essential genes (NEGs), see below.
#' If \code{norm.type.field=NULL}, default edgeR normalization is used instead based on \code{\link{calcNormFactors}}.
#' If \code{norm.type.field=NA}, no normalization is performed beyond library size normalization.
#' \item Consolidation of per-barcode results into per-gene results, in experiments where multiple barcodes (i.e., guides, shRNAs) target the same gene; see below.
#' If \code{gene.field=NA}, no consolidation is performed.
#' }
#'
#' @section Filtering on reference samples:
#' We compute the log-average CPM across reference samples for each barcode.
#' We define a filtering threshold based on the median absolute deviation (MAD) below the median of the log-average CPMs across barcodes.
#' Barcodes that have log-average CPMs that are more than three MADs below the median are discarded,
#' as these correspond to barcodes that were most likely missing from the original pool during its manufacture.
#' 
#' Technically, this filtering strategy is not independent as it can theoretically enrich for barcodes that are upregulated in the reference.
#' In practice, this is not an issue as the filtering is only designed to remove outliers.
#' The density of barcodes around the filter boundary is so low that there will be no noticeable effect on type I error.
#'
#' @section Normalization with control barcodes:
#' Restricting the normalization to only use control barcodes weakens the assumption of a non-DE majority of barcodes.
#' We currently use TMM normalization (i.e., robust average of a ratio) rather than methods based on a robust average of expression within each sample.
#' This avoids inflated errors due to the greater variance of the distribution of expression values compared to the ratio.
#'
#' @section Consolidating barcodes to genes:
#' When \code{method="holm-mid"}, we use the \dQuote{minimum-Holm} method in \code{\link{combineBarcodeTests}} to combine p-values from multiple barcodes into a single p-value per gene.
#' This applies a Holm-Bonferroni correction across all barcodes for each gene before taking the \eqn{k}-th largest p-value for that gene.
#' The idea is to detect genes for which at least \eqn{k} barcodes are significantly differentially abundant, requiring some level of agreement between barcodes to reduce the influence of off-target guides.
#' We also report the log-counts-per-million and log-fold changes of the best barcode within each gene.
#'
#' When \code{method="simes"}, we use Simes' method in \code{\link{combineBarcodeTests}} to combine barcode-level p-values.
#' This is fast, statistically rigorous and robust to correlations between guides for the same gene.
#' It is able to detect genes with only a minority of differentially abundant barcodes, though it will favor genes with many significant barcodes.
#' This tends to be the most sensitive approach but is also more susceptible to off-target effects where one outlier guide drives the gene-level result.
#' Again, we report the log-CPM and log-fold changes of the best barcode within each gene.
#' 
#' When \code{method="fry"}, we use the \code{\link{fry}} method to combine per-barcode statistics into a per-gene p-value.
#' This effectively applies gene set testing machinery on the set of guides for each gene.
#' The aim is to identify genes for which the mean t-statistic is significantly non-zero, which favors consistency among guides more than \code{method="simes"},
#' This may yield low p-values for genes with weak but consistent DE in each guide, at the cost of genes that have strong DE but only in a few guides.
#' Note that this mode is not compatible with ANOVA-like contrasts or with non-zero \code{lfc}.
#'
#' @author Aaron Lun
#'
#' @examples
#' # Setting up a new project.
#' library(gp.sa.core)
#' proj <- tempfile()
#' newProject(proj, title="Aaron's project",
#'     description="This is one of Aaron's projects")
#'
#' # Mocking up an example dataset.
#' N <- 1000
#' mu <- 2^runif(N, 2, 10)
#' y <- matrix(rnbinom(N * 6, mu=mu, size=10), ncol=6)
#' g <- gl(2, 3)
#' 
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(y)
#' rowData(se)$type <- sample(c("endog", "NTG", "NEG"), N, replace=TRUE)
#' rowData(se)$gene <- paste0("GENE_", sample(100, N, replace=TRUE))
#' colData(se)$group <- g
#'
#' report <- file.path(proj, "voom-screen.Rmd")
#' out <- runVoomScreen(se, groups="group",
#'     comparisons=list(c("2", "1")),
#'     norm.type.field="type", norm.type.level="NEG",
#'     reference.field="group", reference.level="1",
#'     gene.field="gene", fname=report)
#'
#' # Returns the Rmarkdown report and analysis results.
#' file.exists(report)
#' out
#' out[[1]]
#'
#' @export
#' @importFrom gp.sa.diff .runVoomCore .defaultEdgeRFilter .defaultEdgeRNormalize .createContrasts
#' @importFrom gp.sa.core .reportStart .reportEnd .createTempRmd .knitAndWrite
#' @importFrom grDevices pdf dev.list dev.off
#' @importFrom methods as
#' @importFrom S4Vectors List
#' @importFrom limma eBayes treat
runVoomScreen <- function(se, groups, comparisons, 
    reference.field, reference.level, norm.type.field, norm.type.level, gene.field, method=c("simes", "holm-mid", "fry"),
    ..., annotation=NULL, lfc=0, robust=TRUE, dup.cor=NULL, contrasts.fun=NULL,
    dump.norm=NULL, fname='voom-screen.Rmd', commit="auto", save.all=NULL)
{
    # Disable graphics devices to avoid showing a whole bunch of plots.
    if (is.null(dev.list())) {
        pdf(file=NULL)
        on.exit(dev.off())
    }

    holding <- .createTempRmd(fname)
    on.exit(unlink(holding))

    .reportStart(fname,
        title="Differential abundance analysis of barcode count data with `voom` and _limma_",
        author=list(list(name="Aaron Lun", affiliation="Genentech gRED B&CB", email="luna@gene.com")),
        call=match.call(),
        commit=commit,
        temporary=holding
    )

    # Setting up all the parts of the analysis that are screen-specific.
    if (is.null(reference.field)) {
        filt <- .defaultEdgeRFilter("barcodes")
    } else if (is.na(reference.field)) {
        filt <- ""
    } else {
        filt <- .filterReference(reference.field, reference.level)
    }

    if (is.null(norm.type.field)) {
        norm <- .defaultEdgeRNormalize("barcodes")
    } else if (is.na(norm.type.field)) {
        norm <- ""
    } else {
        norm <- .normalizeControls(norm.type.field, norm.type.level)
    }
    if (!is.null(dump.norm)) {
        norm <- paste0(norm, sprintf("\n\nWe dump out the normalized abundances in a CSV file for diagnostic use elsewhere.
Do **NOT** use this file as input in other GPSA pipelines.

```{r}
write.csv(file=%s, cpm(y, log=TRUE, prior.count=3))
```", deparse(dump.norm)))
    }

    env <- .runVoomCore(holding, se, groups=groups, comparisons=comparisons,
        lfc=lfc, robust=robust, dup.cor=dup.cor, 
        annotation=annotation, ..., contrasts.fun=contrasts.fun,
        filter=filt, normalize=norm, 
        feature=c("barcode", "barcodes"), analysis="abundance",
        diagnostics=.screen_edgeR_diag_plots(norm.type.field, norm.type.level),
        skip.contrasts=TRUE
    )

    contrast.cmds <- .createContrasts(comparisons, contrasts.fun)
    .screen_contrast_prep(holding, env, gene.field, method=match.arg(method), annotation=annotation)
    .screen_contrast_loop(holding, env, gene.field, lfc=lfc, robust=robust, dup.cor=dup.cor, 
        method=match.arg(method), contrast.cmds=contrast.cmds)

    # Deciding whether or not we can save stuff.
    do.genes <- !is.na(gene.field)
    if (is.null(save.all)) {
        save.all <- !is(se, "SummarizedExperiment") || !is.null(trackinfo(se)$origin)
    }
    if (!save.all) {
        saveable <- NULL
    } else {
        if (do.genes) {
            saveable <- c("barcode.results", "gene.results")
        } else {
            saveable <- "barcode.results"
        }
    }

    .reportEnd(fname, msg="Created report with runVoomScreen().", 
        commit=commit, env=env, to.save.list=saveable, temporary=holding)

    output <- List(barcode=env$barcode.results)
    if (do.genes) {
        output$gene <- env$gene.results
    }
    output
}

####################################
####################################

#' @importFrom gp.sa.diff .defaultEdgeRMDS .defaultEdgeRMD
.screen_edgeR_diag_plots <- function(norm.type.field, norm.type.level) {
    if (is.null(norm.type.field)) {
        md.code <- .defaultEdgeRMD("barcodes")
    } else if (is.na(norm.type.field)) {
        md.code <- ""
    } else {
        md.code <- sprintf("We examine the performance of normalization by creating MD (mean-difference) plots.
Most of the negative control barcodes should have a log-fold change of zero between samples if normalization was successful.

```{r, fig.wide=TRUE, fig.asp=ceiling(ncol(y)/3)/3}
norm.use <- rowData(se)[[%s]][y$genes$origin] %%in%% %s
lcpm <- cpm(y, log=TRUE, prior.count=3)
n <- ncol(lcpm)
par(mfrow=c(ceiling(n/3), 3))
for (i in seq_len(n)) {
    plotMD(lcpm, column=i, status=norm.use)
    abline(h=0, col='red', lty=2)
}
```", deparse(norm.type.field), deparse(norm.type.level))
    }

    paste("# Making diagnostic plots",
        md.code,
        .defaultEdgeRMDS("barcodes"),
        sep="\n\n")
}

####################################
####################################

#' @importFrom gp.sa.core .openChunk .closeChunk .justWrite .evalAndWrite
.screen_contrast_prep <- function(fname, env, gene.field, method="simes", annotation=NULL) {
    do.genes <- !is.na(gene.field)
    if (do.genes) {
        gene_txt <- "\nWe also create a function to convert per-barcode results into per-gene results."
    } else {
        gene_txt <- ""
    }

    .justWrite(fname, sprintf('# Setting up contrasts

## Overview

We set up a function to format the results into a table that lists all the barcodes that we previously filtered out with `NA` statistics,
so as to distinguish between filtered barcodes and those that were not in `se` in the first place.%s', gene_txt))

    # Defining formatting functions for the raw limma output.
    .openChunk(fname)

#    if (do.genes) {
#        listing <- union(gene.field, annotation)
#    } else {
        listing <- annotation
#    }
    .evalAndWrite(fname, env, sprintf('library(gp.sa.core)
library(gp.sa.diff)
library(gp.sa.screen)
barcode_formatter <- function(res, ...) {
    res <- cleanDataFrame(res, se, subset=res$origin,
        anno.fields=%s)
    res$origin <- NULL
    res
}', paste(deparse(listing), collapse="\n        ")))

    if (do.genes) {
        .evalAndWrite(fname, env, sprintf("\ngene.ids <- rowData(se)[[%s]]
gene_formatter <- function(gres) {
    m <- match(rownames(gres), gene.ids)
    cbind(gres, rowData(se)[m,%s,drop=FALSE])
}", deparse(gene.field), paste(deparse(annotation), collapse="\n        ")))
    }
    .closeChunk(fname)

    # Setting up output lists.
    .justWrite(fname, "We set up a `List` to hold all of our output results.

```{r}")
    .evalAndWrite(fname, env, "barcode.results <- List()")
    if (do.genes) {
        .evalAndWrite(fname, env, "gene.results <- List()")
    }
    .closeChunk(fname)
}

#' @importFrom gp.sa.core .openChunk .closeChunk .justWrite .knitAndWrite .evalAndWrite
.screen_contrast_loop <- function(fname, env, gene.field, method="simes", lfc=0, robust=TRUE, dup.cor=NULL, contrast.cmds) 
# Code mostly lifted from limma_contrast_chunk, which was necessary
# to smoothly handle the barcode -> gene result conversion.
{
    do.genes <- !is.na(gene.field)
    extra_eb_code <- if (robust) ", robust=TRUE" else ""
    shrink_eb_cmd <- sprintf("eBayes(fit2%s)", extra_eb_code)
    fry_params <- .get_fry_params(robust=robust, dup.cor=dup.cor)

    for (con in seq_len(nrow(contrast.cmds))) {
        vname <- contrast.cmds$name[con]
        .justWrite(fname, sprintf('## %s', vname))
        .justWrite(fname, "### Setting up the contrast")

        if (!is.null(env$con)) { # clearing out existing 'con'.
            rm("con", envir=env)
        }

        .openChunk(fname)
        .evalAndWrite(fname, env, contrast.cmds$contrast[con])

        if (is.null(env$con)) {
            stop("contrast commands should define a 'con' object")
        }

        gp.sa.diff:::.enforce_contrast_colnames(fname, env)
        .justWrite(fname, "con", trail=FALSE)
        .closeChunk(fname)

        solo <- is.null(dim(env$con)) || ncol(env$con)==1L
        if (!solo) {
            null <- 0
        } else {
            null <- lfc
        }
        
        if (null==0) {
            .justWrite(fname, "### Testing for any differenc

We test the filtered barcodes for significant differences in abundance.")
            cur_shrink_cmd <- shrink_eb_cmd
        } else {
            dlfc <- deparse(null)
            .justWrite(fname, sprintf("### Testing for log-FCs above %s

We test the filtered barcodes for log-fold changes that are significantly more extreme than %s using the `treat` function", dlfc, dlfc))
            cur_shrink_cmd <- sprintf("treat(fit2, lfc=%s%s)", dlfc, extra_eb_code)
        }

        .knitAndWrite(fname, env, sprintf('```{r}
fit2 <- contrasts.fit(fit, con)
fit2 <- %s
res <- topTable(fit2, n=Inf, sort.by="none")
res <- barcode_formatter(res)
head(res[order(res$PValue),])
```', cur_shrink_cmd))

        sum.cmd <- if (solo) "summary(decideTests(fit2))" else "summary(res$FDR <= 0.05)"
        .knitAndWrite(fname, env, sprintf('We report summary statistics for this comparison, defining significant differences in abundance at a FDR threshold of 5%%.

```{r}
%s
```', sum.cmd))

        .knitAndWrite(fname, env, sprintf("We save the results in our output `List` for later use.

```{r}
con.desc <- %s
barcode.results[[con.desc]] <- DAScreenStatFrame(res, se, contrast=con, 
    description=con.desc, method='voom', feature='barcode')
```", deparse(vname)))

        if (do.genes) {
            if (method=="simes") {
                .knitAndWrite(fname, env, "We also consolidate per-barcode statistics into per-gene results using Simes' method.
This tests the global null hypothesis that all barcodes for a gene are not differentially abundant.
We also add the statistics for the best barcode (i.e., that with the lowest $p$-value) for reporting purposes.

```{r}
gres <- combineBarcodeTests(res, gene.ids)
gres <- gene_formatter(gres)
head(gres[order(gres$PValue),])
```")
            } else if (method=="holm-mid") { 
                .knitAndWrite(fname, env, "We also consolidate per-barcode statistics into per-gene results using \"minimum Holm\" method.
This tests the global null hypothesis that fewer than $X$ barcodes for a gene are differentially abundant, where $X$ is defined by `min.sig.n` and `min.sig.prop`.
We also add the statistics for the best barcode (i.e., that with the lowest $p$-value) for reporting purposes.

```{r}
gres <- combineBarcodeTests(res, gene.ids, method='holm-min',
    min.sig.n=3, min.sig.prop=0.4)                                
gres <- gene_formatter(gres)
head(gres[order(gres$PValue),])
```")
            } else {
                .knitAndWrite(fname, env, sprintf("We also consolidate per-barcode statistics into per-gene results using the `fry` method.
This uses the machinery for performing a self-contained gene set test, and applies it to 'barcode sets' corresponding to a single gene.
The null hypothesis is that no barcodes for a given gene are differentially abundant, but favoring genes where most barcodes are DA.
We also add the median log-FC and log-CPM across all barcodes for the each gene.

```{r}
gres <- barcodeSetTest(v, gene.ids, design=design, contrast=con, 
    subset=v$genes$origin, stats=res%s)
gres <- gene_formatter(gres)
head(gres[order(gres$PValue),])
```", fry_params))
            }

            sum.cmd <- if (solo) "table(Sig=gres$FDR <= 0.05, Direction=gres$Direction)" else "summary(gres$FDR <= 0.05)"
            .knitAndWrite(fname, env, sprintf('We report summary statistics for this comparison at the gene level.

```{r}
%s
```', sum.cmd))

            .knitAndWrite(fname, env, "We then save it into our result `List`.

```{r}
gene.results[[con.desc]] <- DAScreenStatFrame(gres, se, contrast=con, 
    description=con.desc, method='voom', feature='gene')
```")
        }
    }

    invisible(NULL)
}

.get_fry_params <- function(dup.cor, robust) {
    included <- character(0)
    if (!is.null(dup.cor)) {
        included <- c(included, "block=dc.block", "correlation=dc$consensus")
    }
    if (robust) {
        included <- c(included, "robust=TRUE")
    }
    if (length(included)) {
        arg.string <- paste(included, collapse=", ")
        arg.string <- paste0(",\n    ", arg.string)
    } else {
        arg.string <- ""
    }
    arg.string
}
