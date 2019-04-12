#' Combine combinatorial barcode counts
#'
#' Combine counts for combinatorial barcodes from multiple files into a single count matrix.
#'
#' @param ... Any number of \linkS4class{DataFrame}s produced by \code{\link{countComboBarcodes}}.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{combination}, a DataFrame containing all unique combinatorial barcodes observed in any \code{...}.
#' Each row corresponds to a barcode and each column contains an integer identifier for the sequence in the variable region.
#' \item \code{counts}, a matrix with number of columns equal to number of objects in \code{...}. 
#' Each row corresponds to a unique combinatorial barcode in \code{keys} and each column represents the count of that barcode in each file of \code{files}.
#' }
#' 
#' @author Aaron Lun
#' @examples
#' df1 <- DataFrame(combination=I(DataFrame(X=1:4, Y=1:4)),
#'    count=sample(10, 4))
#' 
#' df2 <- DataFrame(combination=I(DataFrame(X=1:4, Y=4:1)),
#'    count=sample(10, 4))
#'
#' df3 <- DataFrame(combination=I(DataFrame(X=1, Y=1)),
#'    count=sample(10, 1))
#'
#' combineComboCounts(df1, df2, df3)
#' 
#' @export
combineComboCounts <- function(...) {
    everything <- list(...)
    combined <- do.call(rbind, everything)
    combined$origin <- rep(seq_along(everything), vapply(everything, nrow, FUN.VALUE=0L))

    o <- do.call(order, as.list(combined$combination))
    combined <- combined[o,]
    if (nrow(combined)) {
        any.diff <- lapply(as.list(combined$combination), FUN=function(x) c(TRUE, diff(x)!=0L))
        is.unique <- Reduce("|", any.diff)
    } else {
        is.unique <- logical(0)
    }

    out.counts <- matrix(0L, sum(is.unique), length(everything))
    id <- cumsum(is.unique)
    out.counts[cbind(id, combined$origin)] <- combined$count
    list(combination=combined$combination[is.unique,], counts=out.counts)
}
