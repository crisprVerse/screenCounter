#' Filter on the reference
#'
#' Generate code to compute a filter threshold for removing low-abundance barcodes,
#' based on the average log-abundance of all barcodes in the reference samples.
#'
#' @param ref.field String containing the name of the column of \code{colData} specifying the type of each sample.
#' @param to.use Character vector containing the types of samples that are references.
#'
#' @return 
#' A string containing commands to filter barcodes based on their abundances in the reference samples.
#'
#' @details
#' The output commands assume that there is a \linkS4class{SummarizedExperiment} object named \code{se} and a DGEList object named \code{y} in the evaluation environment.
#' Filtering will be applied based on the reference samples to modify \code{y} in place.
#' 
#' @author Aaron Lun
#' @examples
#' cat(filterReference("condition", "ref"))
#' cat(filterReference("time", 0))
#'
#' @seealso
#' \code{\link{runVoomScreen}}, in which this function is called.
#'
#' @export
#' @rdname filterReference
#' @importFrom edgeR aveLogCPM
.filterReference <- function(ref.field, to.use) {
    sprintf("We use the reference samples from the screen to remove low-abundance barcodes.
These usually represent barcodes that were missing from the original pool (e.g., due to defects in manufacture) and are not of interest.
We define a threshold on the log-average CPM by defining small outliers based on the median absolute deviation.

```{r, fig.cap='Distribution of log-abundances across all barcodes. The blue line represents the median while the red line represents the chosen threshold.'}
is.ref <- colData(se)[[%s]] %%in%% %s
ref.ab <- aveLogCPM(y[,is.ref])

med.ab <- median(ref.ab)
med.ab
mad.ab <- mad(ref.ab, center=med.ab)
threshold <- med.ab - 3*mad.ab
threshold

hist(ref.ab, col='grey80', breaks=20,
    xlab='Average abundance')
abline(v=med.ab, col='blue', lty=2)
abline(v=threshold, col='red', lty=2)
```

It is then straightforward to apply the filter thresholds to the `DGEList` object.

```{r}
filtered <- ref.ab >= threshold
y <- y[filtered,]
summary(filtered)
```", deparse(ref.field), deparse(to.use))
}
