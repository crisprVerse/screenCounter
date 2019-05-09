#' Use \code{voom} to analyze a screen
#'
#' Perform a \code{voom} analysis on a count matrix from a sequencing screen to detect differentially abundant barcodes across samples.
#'
#' @param se A \linkS4class{SummarizedExperiment} object.
#' Alternatively, a string containing a file path to a serialized instance of a SummarizedExperiment.
#' @param ... Further arguments to pass to \code{\link{runVoom}}. 
#' @param reference.field String specifying the column of \code{colData(se)} containing the type of each sample (i.e., reference or not).
#' @param reference.level Character vector specifying the reference levels of the column named by \code{reference.field}.
#' @param norm.type.field String specifying the field of \code{rowData(se)} containing the gene type for each barcode.
#' @param norm.type.level Character vector specifying the gene types on which to perform normalization.
#' @param gene.field String specifying the field of \code{rowData(se)} that contains the gene identifier for each barcode.
#' @param fname String containing the path to an output Rmarkdown file.
#' @param res.dir String containing the name of the directory for result files.
#' @param obj.dir String containing the name of the directory for object files.
#'
#' @return A list containing \code{objects} and \code{results}.
#' \code{objects} is a list containing the \code{fit} object, the MArrayLM object produced by \code{eBayes}.
#' \code{results} is a named list containing \linkS4class{DataFrame}s of result tables from all contrasts.
#'
#' @details
#' This function is largely a wrapper around \code{\link{runVoom}}, with three additional features:
#' \itemize{
#' \item Filtering based on reference samples; see \code{\link{filterReference}}.
#' If \code{reference.field=NULL}, default edgeR filtering is used instead; see \code{\link{defaultEdgeRFilter}}.
#' If \code{reference.field=NA}, no filtering is performed.
#' \item Normalization based on non-targeting genes (NTGs) and/or non-essential genes (NEGs); see \code{\link{normalizeControls}}.
#' If \code{norm.type.field=NULL}, default edgeR normalization is used instead; see \code{\link{defaultEdgeRNormalize}}.
#' If \code{norm.type.field=NA}, no normalization is performed beyond library size normalization.
#' \item Consolidation of per-barcode results into per-gene results, in experiments where multiple barcodes (i.e., guides, shRNAs) target the same gene; see \code{\link{consolidateGenes}}.
#' If \code{gene.field=NA}, no consolidation is performed.
#' }
#'
#' @author Aaron Lun
#'
#' @examples
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
#' tmp <- tempfile(fileext=".Rmd")
#' out <- runVoomScreen(se, groups="group",
#'     comparisons=list(c("2", "1")),
#'     norm.type.field="type", norm.type.level="NEG",
#'     reference.field="group", reference.level="1",
#'     gene.field="gene", fname=tmp)
#' file.exists(tmp)
#'
#' head(out$results$`2-1:de:barcode`)
#' head(out$results$`2-1:de:gene`)
#'
#' @export
#' @importFrom gp.sa.diff runVoomCore createDesignMatrix createContrasts
#' defaultEdgeRFilter defaultEdgeRNormalize
#' @importFrom gp.sa.core makeFrontMatter knitAndWrite newReportPath 
#' @importFrom grDevices pdf dev.list dev.off
runVoomScreen <- function(se, ..., 
    reference.field, reference.level, norm.type.field, norm.type.level, gene.field,
    fname=NULL, res.dir="results", obj.dir="objects")
{
    # Disable graphics devices to avoid showing a whole bunch of plots.
    if (is.null(dev.list())) {
        pdf(file=NULL)
        on.exit(dev.off())
    }

    # Do this before changing wd.
    if (is.character(se)) se <- normalizePath(se) 

    # Choosing whether to output just barcode results, or to consolidate.
    if (is.na(gene.field)) {
        postcon <- NULL
        res.fmt <- "%s"
    } else {
        postcon <- consolidateGenes(gene.field, res.dir)
        res.fmt <- "%s:barcode"
    }

    # Define the output file and temporarily change directory to it.
    # Note that CD'ing is done *after* defining input paths.
    fname <- newReportPath(fname, "voom-screen")
    old <- getwd()
    setwd(dirname(fname))
    fname <- basename(fname)
    on.exit(setwd(old))

    # Using write() here to force overwrite!
    write(file=fname, c(makeFrontMatter(
        title="Differential abundance analysis of barcode count data with `voom` and _limma_",
        author=list(list(name="Aaron Lun", affiliation="Genentech gRED B&CB", email="luna@gene.com"))
    ), ""))

    # Setting up all the parts of the analysis that are screen-specific.
    if (is.null(reference.field)) {
        filt <- defaultEdgeRFilter("barcodes")
    } else if (is.na(reference.field)) {
        filt <- ""
    } else {
        filt <- filterReference(reference.field, reference.level)
    }

    if (is.null(norm.type.field)) {
        norm <- defaultEdgeRNormalize("barcodes")
    } else if (is.na(norm.type.field)) {
        norm <- ""
    } else {
        norm <- normalizeControls(norm.type.field, norm.type.level)
    }

    env <- runVoomCore(fname, se, ..., 
        filter=filt, normalize=norm, 
        feature=c("barcode", "barcodes"), analysis="abundance",
        diagnostics=.screen_edgeR_diag_plots(norm.type.field, norm.type.level),
        post.contrast=postcon,
        res.dir=res.dir, obj.dir=obj.dir
    )

    # Cleaning up.
    y.name <- file.path(obj.dir, "y.rds")
    fit.name <- file.path(obj.dir, "fit.rds")
    knitAndWrite(fname, env, sprintf("# Wrapping up

We save the `y` and `fit` objects to file for later use.

```{r}
saveRDS(y, file=%s)
saveRDS(fit, file=%s)
```

<!-- GPSA_OUTPUT
- path: %s
- path: %s
-->

We also report the session information for our records.

```{r}
sessionInfo()
```", deparse(y.name), deparse(fit.name), y.name, fit.name))

    # Reporting the results.
    all.res <- .retrieve_saved_results(fname)
    res <- lapply(all.res, readRDS)
    names(res) <- sub("\\.rds$", "", basename(all.res))
    list(objects=env, results=res)
}

#' @importFrom gp.sa.diff defaultEdgeRMDS defaultEdgeRMD
.screen_edgeR_diag_plots <- function(norm.type.field, norm.type.level) {
    if (is.null(norm.type.field)) {
        md.code <- defaultEdgeRMD("barcodes")
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
        defaultEdgeRMDS("barcodes"),
        sep="\n\n")
}

# Copied verbatim from gp.sa.diff, couldn't be bothered to export this...
.retrieve_saved_results <- function(fname, pattern="^saveRDS\\((res|stats), file=\"(.*)\"\\)$") {
    read.back <- readLines(fname)
    read.back <- read.back[grep(pattern, read.back)]
    sub(pattern, "\\2", read.back)
}
