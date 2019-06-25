#' DataFrame for differential abundance results
#'
#' The DBAStatFrame and DGAStatFrame classes are literally identical to a standard \linkS4class{DataFrame},
#' and can be used as such in all applications.
#' It is intended to hold results from a differential gene/barcode abundance analyses of sequencing screen data.
#' 
#' @section Constructor:
#' \code{DBAStatFrame(x, parent, ...)} will return a DBAStatFrame object, given the arguments:
#' \itemize{
#' \item \code{x}, a \linkS4class{DataFrame} object or something that can be coerced into one.
#' \item \code{parent}, a \linkS4class{Vector}-like object from which \code{x} was derived.
#' This is usually a \linkS4class{SummarizedExperiment}. 
#' \item \code{...}, other named fields to add to the provenance tracking information via \code{\link{trackinfo}}.
#' }
#'
#' Provenance tracking information is added from \code{trackinfo(parent)} to the output DBAStatFrame.
#' Any additional fields from \code{...} are also added, overwriting existing fields if they have the same name.
#'
#' The \code{DGAStatFrame} constructor has the same behavior, only differing in that it returns a DGAStatFrame. 
#' The DGAStatFrame and DBAStatFrame subclasses are functionally equivalent but the former is intended to hold gene-level results while the latter is intended to hold barcode-level results.
#'
#' @section Checking metadata:
#' \code{trackcheck(x)} will check for the presence of correct provenance fields in a DBAStatFrame or DGAStatFrame \code{x}
#' and return a list of this information.
#'
#' In addition to the provenance fields required by \linkS4class{DiffStatFrame}, we also require
#' \code{method}, a string specifying the differential abundance analysis method that was used (currently only \code{"voom"}).
#'
#' A \code{type} field specifying the result type in a controlled vocabulary is included in the output,
#' and will overwrite any \code{type} field in \code{trackinfo(x)}.
#'
#' @author Aaron Lun
#' @examples
#' library(gp.sa.core)
#'
#' # Mocking up an input into runVoomScreen()
#' library(SummarizedExperiment)
#' se.input <- SummarizedExperiment()
#' trackinfo(se.input)$origin <- list(id="SOME_ID")
#'
#' # Mocking up the aftermath of a DE analysis:
#' de.output <- DataFrame(N=1:10, PValue=0:9/10)
#'
#' Y <- DBAStatFrame(de.output, se.input,
#'     contrast=c(A=1, B=-1),
#'     method="voom", description="I did voom")
#' Y 
#' 
#' @seealso
#' \linkS4class{DiffStatFrame}, from which this class is derived.
#' 
#' @docType class
#' @name DBAStatFrame
#' @aliases DBAStatFrame DBAStatFrame-class trackcheck,DBAStatFrame-method
#' DGAStatFrame DGAStatFrame-class trackcheck,DGAStatFrame-method
NULL

#' @export
#' @importFrom gp.sa.core createTSFSubclass
DBAStatFrame <- function(x, parent, ...) {
    x <- createTSFSubclass(x, parent, ...)
    as(x, "DBAStatFrame")
}

#' @export
#' @importFrom gp.sa.core createTSFSubclass
DGAStatFrame <- function(x, parent, ...) {
    x <- createTSFSubclass(x, parent, ...)
    as(x, "DGAStatFrame")
}

#' @export
#' @importFrom gp.sa.core quickError trackcheck
setMethod("trackcheck", "DBAStatFrame", function(x) {
    out <- callNextMethod()
    .common_check(x, out)
    out$type <- "Differential barcode abundance result"
    out
})

#' @export
#' @importFrom gp.sa.core quickError trackcheck
setMethod("trackcheck", "DGAStatFrame", function(x) {
    out <- callNextMethod()
    .common_check(x, out)
    out$type <- "Differential gene abundance result"
    out
})

.common_check <- function(x, info) {
    if (!"method" %in% names(info)) {
        quickError(x, "method", "a string specifying the differential analysis method that was used")
    }
    if (!info$method %in% c("voom")) {
        stop("'method' should be 'voom'")
    }
    invisible(NULL)
}
