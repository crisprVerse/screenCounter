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
#' Normalization factors will be computed using a median-based approach and assigned back to \code{y} for downstream analysis.
#'
#' Here, we use a custom normalization approach that is, unlike conventional methods for RNA-seq data, \emph{not} based on ratios.
#' Rather, the median is directly computed from the counts of the (control) barcodes in each sample.
#' This is a simpler approach that can be more precise than the use of ratios, under the assumption that the barcodes specified by \code{to.use} are present in equimolar amounts.
#' 
#' (In the presence of a large proportion of differential features, the error of the median will be proportional to the variance of the distribution.
#' Equal molarity means that the variance of the counts is smaller than that of the ratios, as the latter effectively involves adding the variance of the log-counts.
#' On the other hand, our approach will be more inaccurate than ratios when the abundance distribution has large variance, e.g., for bulk RNA-seq data.)
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
#' @importFrom edgeR calcNormFactors
normalizeControls <- function(type.field, to.use) {
    normalize <- sprintf("We define a normalization factor for each sample by taking the median abundance of a set of control barcodes.
This avoids composition biases due to genuine changes in abundance for barcodes associated with relevant biological processes.

```{r}
norm.use <- rowData(se)[[%s]][y$origin] %%in%% %s
to.use <- y$counts[norm.use,,drop=FALSE]
scaling <- apply(to.use, 2, median)
new.nf <- scaling / y$samples$lib.size
y$samples$norm.factors <- new.nf/mean(new.nf)", deparse(type.field), deparse(to.use))

    normalize <- c(normalize, "head(y$samples, 10)
summary(y$samples$norm.factors)
```")
    paste(normalize, collapse="\n")
}
