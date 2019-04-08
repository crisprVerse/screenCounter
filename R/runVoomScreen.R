#' Use \code{voom} to analyze a screen
#'
#' Perform a \code{voom} analysis on a count matrix from a sequencing screen to detect differentially abundant barcodes across samples.
#'
#' @param se A \linkS4class{SummarizedExperiment} object.
#' Alternatively, a string containing a file path to a serialized instance of a SummarizedExperiment.
#' @param groups String specifying the column of \code{colData(se)} containing the grouping factor, see \code{\link{createDesignMatrix}}.
#' @param comparisons A list of character vectors specifying the comparisons to perform between groups, see \code{\link{createContrasts}}.
#' @param covariates Character vector specifying the columns of \code{colData(se)} containing continuous covariates of interest,  see \code{\link{createDesignMatrix}}.
#' @param block Character vector specifying additional blocking factors or covariates that are \emph{not} of interest, see \code{\link{createDesignMatrix}}.
#' @param ... Further arguments to pass to \code{\link{runVoom}}. 
#' @param type.field String specifying the field of \code{rowData} specifying the gene type for each barcode.
#' @param norm.types Character vector specifying the gene types on which to perform normalization.
#' @param design.fun A function to create a custom design matrix, see \code{?\link{createDesignMatrix}}.
#' @param contrasts.fun A list of custom contrasts information, including contrast-generating functions; see \code{?\link{createContrasts}}.
#' @param fname String containing the path to an output Rmarkdown file.
#'
#' @return A list containing \code{objects} and \code{results}.
#' \code{objects} is a list containing the \code{fit} object, the MArrayLM object produced by \code{eBayes}.
#' \code{results} is a named list containing \linkS4class{DataFrame}s of result tables from all contrasts.
#'
#' @details
#' This function is largely a wrapper around \code{\link{runVoom}}, with two additional features:
#' \itemize{
#' \item Normalization based on non-targeting genes (NTGs) and/or non-essential genes (NEGs).
#' These are negative controls that may yield more accurate normalization factors when the majority of other genes are DE.
#' \item Consolidation of per-barcode results into per-gene results, in experiments where multiple barcodes (i.e., guides, shRNAs) target the same gene.
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
#' colData(se)$group <- g
#'
#' out <- runVoomScreen(se, groups="group",
#'     comparisons=list(c("2", "1")),
#'     type.field="type", norm.types="NEG")
#' head(out$results$`2-1`)
#'
#' @export
#' @importFrom gp.sa.diff runVoomCore createDesignMatrix createContrasts
#' @importFrom gp.sa.core makeFrontMatter pathFromRoot knitAndWrite newDirectoryPath
runVoomScreen <- function(se, groups, comparisons, covariates=NULL, block=NULL, 
    ..., type.field="XXX", norm.types=NULL,
    design.fun=NULL, contrasts.fun=NULL, fname=NULL)
{
    # Disable graphics devices to avoid showing a whole bunch of plots.
    if (is.null(dev.list())) {
        pdf(file=NULL)
        on.exit(dev.off())
    }

    # Defining inputs and outputs.
    if (is.character(se)) {
        se <- pathFromRoot(se)
    }
    contrast.cmds <- createContrasts(comparisons, contrasts.fun)
    by.barcodes <- paste0(contrast.cmds$name, "_barcode")
    by.genes <- paste0(contrast.cmds$name, "_gene")
    all.results <- c(by.barcodes, by.genes)
    contrast.cmds$name <- by.barcodes

    # Define the output file and temporarily change directory to it.
    # Note that CD'ing is done *after* defining input paths.
    if (is.null(fname)) {
        fname <- file.path(newDirectoryPath("voom-screen"), "report.Rmd")
    }
    old <- getwd()
    setwd(dirname(fname))
    fname <- basename(fname)
    on.exit(setwd(old))

    # Using write() here to force overwrite!
    write(file=fname,
        c(makeFrontMatter(
            title="Differential abundance analysis of barcode count data with `voom` and _limma_",
            author=list(list(name="Aaron Lun", affiliation="Genentech gRED B&CB", email="luna@gene.com")),
            dependencies=if (is.character(se)) se else NULL,
            generated=c(file.path("results", outer(all.results, c(".csv", ".rds"), paste0)), "fit.rds")
        ), "")
    )

    design.cmds <- createDesignMatrix(groups, covariates, block, design.fun)
    design.cmds <- sprintf('## Creating the design matrix

We construct a design matrix based on the sample-specific metadata in our `se` object.
This describes the experimental design to be used for the analysis.

```{r}
%s
colnames(design)
```', design.cmds)

    env <- runVoomCore(se, design.cmds, contrast.cmds, ..., fname=fname,
        filter=gp.sa.diff:::.default_edgeR_filter("genes"), 
        normalize=tmmOnSubset(type.field=type.field, to.use=norm.types),
        out.format=file.path("results", "%s"))

    # Cleaning up.
    knitAndWrite(fname, env, "# Wrapping up

We save the `fit` object to file for later use.

```{r}
saveRDS(fit, file='fit.rds')
```

We also report the session information for our records.

```{r}
sessionInfo()
```")

    # Reporting the results.
    fit <- env$fit
    base.dir <- dirname(fname)
    res <- lapply(all.results, function(x) readRDS(file.path("results", paste0(x, ".rds"))))
    names(res) <- all.results
    list(objects=list(fit=fit), results=res)
}
