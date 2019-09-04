#' DataFrame for differential abundance results
#'
#' The DAScreenStatFrame class is literally identical to a standard \linkS4class{DataFrame},
#' and can be used as such in all applications.
#' It is intended to hold results from a differential abundance analysis of barcode sequencing data from high-throughput CRISPR/siRNA screens.
#' 
#' @section Constructor:
#' \code{DAScreenStatFrame(x, parent, ...)} will return a DAScreenStatFrame object, given the arguments:
#' \itemize{
#' \item \code{x}, a \linkS4class{DataFrame} object or something that can be coerced into one.
#' \item \code{parent}, a \linkS4class{Vector}-like object from which \code{x} was derived.
#' This is usually a \linkS4class{SummarizedExperiment}. 
#' \item \code{...}, other named fields to add to the provenance tracking information via \code{\link{trackinfo}}.
#' }
#'
#' Provenance tracking information is added from \code{trackinfo(parent)} to the output DAScreenStatFrame.
#' Any additional fields from \code{...} are also added, overwriting existing fields if they have the same name.
#'
#' @section Checking metadata:
#' \code{\link{.trackCheck}(x)} will check for the presence of correct provenance fields in a DAScreenStatFrame \code{x}
#' and return a list of this information.
#' Required fields include:
#' \itemize{
#' \item \code{origin}, the sources from which \code{x} was derived - see \linkS4class{TrackedDataFrame}.
#' \item \code{type}, the type of result - see \linkS4class{TrackedDataFrame}.
#' \item \code{description}, a plain-English description - see \linkS4class{TrackedDataFrame}.
#' \item \code{contrast}, the contrast vector/matrix - see \linkS4class{DiffStatFrameBasic}.
#' \item \code{method}, a string specifying the DE analysis method that was used - see \linkS4class{DGEStatFrame}.
#' \item \code{feature}, a string specifying the type of feature in each row of \code{x}.
#' Currently, this can be \code{"gene"} or \code{"barcode"}.
#' }
#' A \code{type} field specifying the \dQuote{differential screen abundance} result type in a controlled vocabulary is included in the output,
#' and will overwrite any \code{type} field in \code{trackinfo(x)}.
#'
#' @section Checking column names:
#' \code{\link{.trackCheck}(x)} will check that the column names of \code{x} include:
#' \itemize{
#' \item \code{"PValue"}, the p-value for each gene.
#' \item \code{"FDR"}, the Benjamini-Hochberg adjusted p-value for each gene.
#' \item \code{"LogCPM"}, the average log2-counts-per-million for each gene.
#' }
#' It should also contain either one \code{"LogFC"} field or multiple \code{"LogFC."}-prefixed fields.
#' These contain the log2-fold changes for a single contrast vector or an ANOVA-like contrast matrix, respectively.
#'
#' @author Aaron Lun
#' @examples
#' library(gp.sa.core)
#'
#' # Mocking up an input into runVoomScreen()
#' library(SummarizedExperiment)
#' se.input <- SummarizedExperiment()
#' trackinfo(se.input)$origin <- list(list(id="SOME_ID"))
#'
#' # Mocking up the aftermath of a DE analysis:
#' de.output <- DataFrame(LogFC=1:10, LogCPM=1:10, 
#'    PValue=0:9/10, FDR=0:9/10)
#'
#' Y <- DAScreenStatFrame(de.output, se.input,
#'     contrast=c(A=1, B=-1), feature="barcode",
#'     method="voom", description="I did voom")
#' Y 
#' 
#' @seealso
#' \linkS4class{DGEStatFrame}, from which this class is derived.
#' 
#' @docType class
#' @name DAScreenStatFrame
#' @aliases DAScreenStatFrame DAScreenStatFrame-class .trackCheck,DAScreenStatFrame-method
NULL

#' @export
#' @importFrom gp.sa.core .createTDFSubclass
DAScreenStatFrame <- function(x, parent, ...) {
    x <- .createTDFSubclass(x, parent, ...)
    as(x, "DAScreenStatFrame")
}

#' @export
#' @importFrom gp.sa.core .trackCheck trackinfo trackinfo<- .quickError
setMethod(".trackCheck", "DAScreenStatFrame", function(x) {
    out <- callNextMethod()
    out$type <- "differential screen abundance"
    if (is.null(out$feature) || !is.character(out$feature) || length(out$feature)!=1L) {
        .quickError(x, "feature", "a string specifying the type of feature in the rows")
    }
    if (!out$feature %in% c("gene", "barcode")) {
        stop("'feature' should be either 'gene' or 'barcode'")
    }
    out
})
