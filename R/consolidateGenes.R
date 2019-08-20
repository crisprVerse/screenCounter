#' Consolidate to genes
#'
#' Consolidate barcode-level statistics to gene-level statistics.
#'
#' @param gene.field String specifying the field of \code{rowData(se)} that contains the gene identifier for each barcode.
#'
#' @return
#' A function that takes a machine-readable name of the contrast (see \code{\link{.createContrasts}}),
#' and returns a Rmarkdown chunk containing code to store per-barcode and per-gene results in the \code{all.results} list.
#'
#' @details
#' The returned function needs a string equivalent to the \code{name} field in the DataFrame from \code{\link{.createContrasts}}.
#' This is used to determine the name of the output files to which to save results.
#'
#' The Rmarkdown chunk expects the SummarizedExperiment \code{se}, a \code{DGEList} object \code{y} 
#' and a per-barcode result table \code{res} to be present in the evaluation environment.
#'
#' Note that all results for a single contrast are saved as a single entry of \code{all.results};
#' per-barcode and per-gene results are accommodated by making this entry a list of length two to hold the separate result tables.
#' It is useful to unroll this list for easier parsing by users later, see \code{\link{runVoomScreen}}.
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
        per.barcode <- sprintf("We save the results in our output `List` for later use.

```{r}                               
con.desc <- %s
res <- DBAStatFrame(res, se, contrast=con, 
    description=con.desc, method=%s)
```", deparse(vname), deparse(method))

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
gres <- DGAStatFrame(gene.stats, se, contrast=con,
    description=con.desc, method=%s)
all.results[[con.desc]] <- list(barcode=res, gene=gres)
```", to.add, deparse(method))

        paste0(per.barcode, "\n\n", per.gene)
    }
}

.default_postcon <- function(vname, method="voom") {
    sprintf("We save the results in our output `List` for later use.

```{r}
con.desc <- %s
res <- DBAStatFrame(res, se, contrast=con, 
    description=con.desc, method=%s)
all.results[[con.desc]] <- res
```", deparse(vname), deparse(method))
}

