#' Consolidate to genes
#'
#' Consolidate barcode-level statistics to gene-level statistics.
#'
#' @param gene.field String specifying the field of \code{rowData(se)} that contains the gene identifier for each barcode.
#'
#' @return
#' A function that takes a machine-readable name of the contrast (see \code{\link{.createContrasts}}),
#' and returns a Rmarkdown chunk containing code to:
#' \itemize{
#' \item store per-barcode results in the \code{all.results} list, named with a \code{:barcode} suffix.
#' \item combine barcode-level statistics for that contrast into gene-level statistics.
#' \item store per-gene results in the \code{all.results} list, named with a \code{:gene} suffix.
#' }
#'
#' @details
#' The returned function needs a string equivalent to the \code{name} field in the DataFrame from \code{\link{.createContrasts}}.
#' This is used to determine the name of the output files to which to save results.
#'
#' The Rmarkdown chunk expects the SummarizedExperiment \code{se}, a \code{DGEList} object \code{y} and a result table \code{res} to be present in the evaluation environment.
#' 
#' @examples
#' FUN <- .consolidateGenes("gene.type")
#' cat(FUN("A-B"))
#' @export
#' @rdname consolidateGenes
#' @importFrom csaw combineTests
.consolidateGenes <- function (gene.field) {
    to.add <- deparse(gene.field)
    function(vname, method="voom") {
        per.barcode <- .default_postcon(vname)
        per.gene <- sprintf("We also consolidate per-barcode statistics into per-gene results using Simes' method.
This tests the joint null hypothesis that all barcodes for a gene are not differentially abundant.
We also add the statistics for the best barcode (i.e., that with the lowest $p$-value) for reporting purposes.

```{r}
gene.ids <- rowData(se)[[%s]]
gene.stats <- barcodes2genes(res, gene.ids)
gene.stats
```

We expand the result table so that there is one row per gene, and save it into our result `List`.

```{r}
con.desc <- %s
all.results[[con.desc]] <- DGAStatFrame(gene.stats, se, contrast=con,
    description=con.desc, method=%s)
```", to.add, deparse(paste(vname, "(gene)")), deparse(method))

        paste0(per.barcode, "\n\n", per.gene)
    }
}

.default_postcon <- function(vname, method="voom") {
    vname <- paste(vname, "(barcode)")
    sprintf("We save the results in our output `List` for later use.

```{r}
con.desc <- %s
res <- DBAStatFrame(res, se, contrast=con, 
    description=con.desc, method=%s)
all.results[[con.desc]] <- res
```", deparse(vname), deparse(method))
}

