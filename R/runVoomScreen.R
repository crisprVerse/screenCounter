#' Use \code{voom} to analyze a screen
#'
#' Perform a \code{voom} analysis on a count matrix from a sequencing screen to detect differentially abundant barcodes across samples.
#'
#' @param se A \linkS4class{SummarizedExperiment} object.
#' Alternatively, a database ID string or list, see \code{?"\link{gp.sa.core-data-inputs}"}.
#' @param ... Further arguments to pass to \code{\link{runVoom}}. 
#' @param reference.field String specifying the column of \code{colData(se)} containing the type of each sample (i.e., reference or not).
#' @param reference.level Character vector specifying the reference levels of the column named by \code{reference.field}.
#' @param norm.type.field String specifying the field of \code{rowData(se)} containing the gene type for each barcode.
#' @param norm.type.level Character vector specifying the gene types on which to perform normalization.
#' @param gene.field String specifying the field of \code{rowData(se)} that contains the gene identifier for each barcode.
#' @param fname String containing the path to an output Rmarkdown file.
#' @param commit String specifying the auto-committing behavior, see \code{?"\link{gp.sa.core-auto-commits}"}.
#' @param save.all Logical scalar indicating whether the returned \linkS4class{DataFrame}s should also be saved to file.
#' Defaults to \code{FALSE} if \code{se} is a SummarizedExperiment, and \code{TRUE} otherwise.
#'
#' @return A \linkS4class{List} containing \linkS4class{DBAStatFrame} and \linkS4class{DGAStatFrame} objects of result tables from all contrasts.
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
#' We use Simes' method as it is fast, statistically rigorous and able to detect genes with only a minority of differentially abundant barcodes.
#' The exact implementation uses the \code{\link{combineTests}} function from the \pkg{csaw} package.
#' We also report the abundance and log-fold changes of the best barcode within each gene.
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
#' @importFrom gp.sa.diff .runVoomCore .defaultEdgeRFilter .defaultEdgeRNormalize .findDFsToSave
#' @importFrom gp.sa.core .reportStart .reportEnd
#' @importFrom grDevices pdf dev.list dev.off
#' @importFrom methods as
runVoomScreen <- function(se, ..., 
    reference.field, reference.level, norm.type.field, norm.type.level, gene.field,
    fname='voom-screen.Rmd', commit="auto", save.all=NULL)
{
    # Disable graphics devices to avoid showing a whole bunch of plots.
    if (is.null(dev.list())) {
        pdf(file=NULL)
        on.exit(dev.off())
    }

    # Choosing whether to output just barcode results, or to consolidate.
    if (is.na(gene.field)) {
        postcon <- .default_postcon
    } else {
        postcon <- .consolidateGenes(gene.field)
    }

    .reportStart(fname,
        title="Differential abundance analysis of barcode count data with `voom` and _limma_",
        author=list(list(name="Aaron Lun", affiliation="Genentech gRED B&CB", email="luna@gene.com")),
        call=match.call(),
        commit=commit
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

    env <- .runVoomCore(fname, se, ..., 
        filter=filt, normalize=norm, 
        feature=c("barcode", "barcodes"), analysis="abundance",
        diagnostics=.screen_edgeR_diag_plots(norm.type.field, norm.type.level),
        post.contrast=postcon
    )

    .reportEnd(fname, msg="Created report with runVoomScreen().", 
        commit=commit, env=env, to.save=.findDFsToSave(se, env, save.all))

    env$all.results
}

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
