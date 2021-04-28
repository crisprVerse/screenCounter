#' Combine combinatorial barcode counts
#'
#' Combine counts for combinatorial barcodes from multiple files into a single count matrix.
#'
#' @param ... Any number of \linkS4class{DataFrame}s produced by \code{\link{countComboBarcodes}}.
#'
#' @return A \linkS4class{DataFrame} containing:
#' \itemize{
#' \item \code{combinations}, a DataFrame containing all unique combinatorial barcodes observed in any \code{...}.
#' Each row corresponds to a barcode and each column contains an identifier (either integer or character) for the sequence in the variable region.
#' \item \code{counts}, a matrix with number of columns equal to number of objects in \code{...}. 
#' Each row corresponds to a unique combinatorial barcode in \code{keys} and each column represents the count of that barcode in each entry if \code{...}.
#' Column names are set to the names of \code{...}, if supplied.
#' }
#' 
#' @author Aaron Lun
#' @examples
#' df1 <- DataFrame(combinations=I(DataFrame(X=1:4, Y=1:4)),
#'    counts=sample(10, 4))
#' 
#' df2 <- DataFrame(combinations=I(DataFrame(X=1:4, Y=4:1)),
#'    counts=sample(10, 4))
#'
#' df3 <- DataFrame(combinations=I(DataFrame(X=1, Y=1)),
#'    counts=sample(10, 1))
#'
#' combineComboCounts(df1, df2, df3)
#' 
#' @export
combineComboCounts <- function(...) {
    everything <- list(...)
    combined <- do.call(rbind, everything)
    combined$origin <- rep(seq_along(everything), vapply(everything, nrow, FUN.VALUE=0L))

    o <- do.call(order, as.list(combined$combinations))
    combined <- combined[o,]
    if (nrow(combined)) {
        any.diff <- lapply(as.list(combined$combinations), FUN=function(x) c(TRUE, x[-1L]!=x[-length(x)]))
        is.unique <- Reduce("|", any.diff)
    } else {
        is.unique <- logical(0)
    }

    out.counts <- matrix(0L, sum(is.unique), length(everything))
    colnames(out.counts) <- names(everything)

    id <- cumsum(is.unique)
    df <- DataFrame(X=id, Y=combined$origin)
    if (anyDuplicated(df)) {
        stop("barcode combinations should be unique within each DataFrame in '...'")
    }

    out.counts[cbind(id, combined$origin)] <- combined$counts

    DataFrame(combinations=I(combined$combinations[is.unique,]), counts=I(out.counts))
}
