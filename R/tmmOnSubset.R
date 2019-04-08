#' Subsetted TMM normalization
#'
#' Generate code to perform TMM normalization on a subset of features,
#' usually barcodes for control genes to satisfy the assumption of a majority of non-DE features. 
#'
#' @param type.field String containing the name of the column of \code{rowData} specifying the type of the gene targeted by each barcode.
#' @param to.use Character vector containing the types to use during normalization.
#' If \code{NULL}, all barcodes are used.
#'
#' @return
#' A string containing commands to compute normalization factors based on the specified subset of barcodes.
#' 
#' @details
#' The output commands assume that there is a \linkS4class{SummarizedExperiment} object named \code{se} and a DGEList object named \code{y} in the evaluation environment.
#' The subsetted normalization factors will be assigned back to \code{y} for downstream analysis.
#'
#' @author Aaron Lun
#' @examples
#' cat(tmmOnSubset("gene.type", c("NTG", "NEG")))
#'
#' @export
#' @importFrom edgeR calcNormFactors
tmmOnSubset <- function(type.field, to.use) {
    normalize <- "We then perform trimmed mean of M-values (TMM) normalization.
This estimates normalization factors for all samples that account for composition biases, i.e., scaling effects beyond library size."

    if (!is.null(to.use)) {
        normalize <- c(normalize, sprintf("Here, we use only a subset of control genes to avoid problems with violations of the assumption of the non-DE majority.
                       
```{r}
norm.use <- rowData(se)[[%s]][filtered] %%in%% %s
new.nf <- calcNormFactors(y[norm.use,])$samples$norm.factors
y$samples$norm.factors <- new.nf", deparse(type.field), deparse(to.use)))
    } else {
        normalize <- c(normalize, "
```{r}
y <- calcNormFactors(y)")
    }

    normalize <- c(normalize, "head(y$samples, 10)
summary(y$samples$norm.factors)
```")
    paste(normalize, collapse="\n")
}
