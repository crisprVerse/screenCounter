#' Consolidate to genes
#'
#' Consolidate barcode-level statistics to gene-level statistics.
#'
#' @param gene.field String specifying the field of \code{rowData(se)} that contains the gene identifier for each barcode.
#'
#' @return
#' A function that takes a row of a \linkS4class{DataFrame} produced by \code{\link{createContrasts}},
#' and returns a Rmarkdown chunk containing code to combine barcode-level statistics for that contrast into gene-level statistics.
#'
#' @details
#' We use Simes' method as it is fast, statistically rigorous and able to detect genes with only a minority of differentially abundant barcodes.
#' The exact implementation uses the \code{\link{combineTests}} function from the \pkg{csaw} package.
#' We also report the abundance and log-fold changes of the best barcode within each gene.
#'
#' The output function technically only needs the \code{name} field in the DataFrame from \code{\link{createContrasts}}.
#' This is used to determine the name of the output files to which to save results.
#'
#' The Rmarkdown chunk expects the SummarizedExperiment \code{se}, a \code{DGEList} object \code{y} and a result table \code{res}.
#' 
#' @examples
#' FUN <- consolidateGenes("gene.type")
#' cat(FUN(DataFrame(name="A-B")))
#' @export
#' @importFrom csaw combineTests
consolidateGenes <- function (gene.field) {
    to.add <- deparse(gene.field)
    function(contrast) {
        vname <- file.path("results", paste0(contrast$name, "_gene"))
        sprintf("### By gene 
                
We also consolidate per-barcode statistics into per-gene results using Simes' method.
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

We expand the result table so that there is one row per gene, and save it to file.

```{r}
stats <- cbind(per.gene, best)
all.genes <- sort(unique(rowData(se)[[%s]]))
expander <- match(all.genes, rownames(stats))
stats <- stats[expander,]
rownames(stats) <- all.genes

write.csv(file=%s, stats)
saveRDS(stats, file=%s)
```", to.add, to.add, deparse(paste0(vname, ".csv")), deparse(paste0(vname, ".rds")))
    }
}
