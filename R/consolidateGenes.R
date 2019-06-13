#' Consolidate to genes
#'
#' Consolidate barcode-level statistics to gene-level statistics.
#'
#' @param gene.field String specifying the field of \code{rowData(se)} that contains the gene identifier for each barcode.
#'
#' @return
#' A function that takes a machine-readable name of the contrast (see \code{\link{createContrasts}}),
#' and returns a Rmarkdown chunk containing code to:
#' \itemize{
#' \item store per-barcode results in the \code{all.results} list, named with a \code{:barcode} suffix.
#' \item combine barcode-level statistics for that contrast into gene-level statistics.
#' \item store per-gene results in the \code{all.results} list, named with a \code{:gene} suffix.
#' }
#'
#' @details
#' We use Simes' method as it is fast, statistically rigorous and able to detect genes with only a minority of differentially abundant barcodes.
#' The exact implementation uses the \code{\link{combineTests}} function from the \pkg{csaw} package.
#' We also report the abundance and log-fold changes of the best barcode within each gene.
#'
#' The returned function needs a string equivalent to the \code{name} field in the DataFrame from \code{\link{createContrasts}}.
#' This is used to determine the name of the output files to which to save results.
#'
#' The Rmarkdown chunk expects the SummarizedExperiment \code{se}, a \code{DGEList} object \code{y} and a result table \code{res} to be present in the evaluation environment.
#' 
#' @examples
#' FUN <- consolidateGenes("gene.type")
#' cat(FUN("A-B"))
#' @export
#' @importFrom csaw combineTests
consolidateGenes <- function (gene.field) {
    to.add <- deparse(gene.field)
    function(vname) {
        per.barcode <- .default_postcon(vname)
        per.gene <- sprintf("We also consolidate per-barcode statistics into per-gene results using Simes' method.
This tests the joint null hypothesis that all barcodes for a gene are not differentially abundant.

```{r}
grouping <- rowData(se)[[%s]][y$genes$origin]
subres <- res[y$genes$origin,]
per.gene <- csaw::combineTests(grouping, subres,
    pval.col='PValue', fc.col=grep('LogFC', colnames(res)))
colnames(per.gene)[1] <- 'nbarcodes'
head(per.gene)
```

We add the statistics for the best barcode (i.e., that with the lowest $p$-value) for reporting purposes.

```{r}
best <- csaw::getBestTest(grouping, subres, pval.col='PValue')
best <- best[,c(which(colnames(best)=='AverageAbundance'), 
    grep('LogFC', colnames(best))),drop=FALSE]
colnames(best) <- paste0('Best', colnames(best))
head(best)
```

We expand the result table so that there is one row per gene.

```{r}
stats <- cbind(per.gene, best)
all.genes <- sort(unique(rowData(se)[[%s]]))
expander <- match(all.genes, rownames(stats))
stats <- stats[expander,]
rownames(stats) <- all.genes
```

... and save it into our result `List`.

```{r}
stats <- as(stats, 'DGEStatFrame')
trackinfo(stats) <- trackinfo(se)
trackinfo(stats)$contrast <- con
con.desc <- %s
trackinfo(stats)$description <- con.desc
all.results[[con.desc]] <- stats
```", to.add, to.add, deparse(paste(vname, "(gene)")))

        paste0(per.barcode, "\n\n", per.gene)
    }
}

.default_postcon <- function(vname) {
    vname <- paste(vname, "(barcode)")
    sprintf("We save the results in our output `List` for later use.

```{r}
res <- as(res, 'DGEStatFrame')
trackinfo(res) <- trackinfo(se)
trackinfo(res)$contrast <- con
con.desc <- %s
trackinfo(res)$description <- con.desc
all.results[[con.desc]] <- res
```", deparse(vname))
}            

