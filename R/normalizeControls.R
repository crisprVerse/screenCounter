#' Normalization on controls
#'
#' Generate code to perform median-based normalization on a subset of features.
#' usually barcodes for control genes to satisfy the assumption of a majority of non-DE features. 
#'
#' @param type.field String containing the name of the column of \code{rowData} specifying the type of the gene targeted by each barcode.
#' @param to.use Character vector containing the types to use during normalization.
#'
#' @return
#' A string containing commands to compute normalization factors based on the specified subset of barcodes.
#' 
#' @details
#' The output commands assume that there is a \linkS4class{SummarizedExperiment} object named \code{se} and a DGEList object named \code{y} in the evaluation environment.
#' Normalization factors are computed using TMM normalization on the subset of control barcodes,
#' and assigned back to \code{y} for downstream analysis.
#'
#' @author Aaron Lun
#' @examples
#' cat(normalizeControls("gene.type", "NTG"))
#' cat(normalizeControls("gene.type", "NEG"))
#'
#' # One can also specify multiple features,
#' # but this is probably unwise.
#' cat(normalizeControls("gene.type", c("NTG", "NEG")))
#'
#' @export
#' @rdname normalizeControls
#' @importFrom edgeR calcNormFactors
.normalizeControls <- function(type.field, to.use) {
    normalize <- sprintf("We define a normalization factor for each sample by taking the median abundance of a set of control barcodes.
This avoids composition biases due to genuine changes in abundance for barcodes associated with relevant biological processes.

```{r}
norm.use <- rowData(se)[[%s]][y$genes$origin] %%in%% %s
ysub <- calcNormFactors(y[norm.use,])
y$samples$norm.factors <- ysub$samples$norm.factors", deparse(type.field), deparse(to.use))

    normalize <- c(normalize, "head(y$samples, 10)
summary(y$samples$norm.factors)
```")
    paste(normalize, collapse="\n")
}
