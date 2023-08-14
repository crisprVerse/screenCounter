#' Match sequences to a pool of barcodes
#'
#' Pretty much what it says on the tin.
#' Useful for matching observed sequences (e.g., from \code{\link{countRandomBarcodes}}) to a pool of known barcode sequences,
#' accounting for substitutions and ambiguous IUPAC codes.
#'
#' @param sequences Character vector of observed sequences.
#' @param choices Character vector of barcode sequences.
#' @param substitutions Integer scalar specifying the maximum number of substitutions when considering a match. 
#' @param reverse Whether to match \code{sequences} to the reverse complement of \code{choices}.
#'
#' @return \linkS4class{DataFrame} with one row per entry of \code{sequences}, containing the following fields:
#' \itemize{
#' \item \code{index}, the index of the matching barcode in \code{choices}.
#' This is set to \code{NA} if no unambiguous match is found.
#' \item \code{mismatches}, the number of mismatching bases with the assigned barcode.
#' This is set to \code{NA} if \code{index} is \code{NA}.
#' }
#'
#' @author Aaron Lun
#'
#' @examples
#' choices <- c("AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT")
#' matchBarcodes(c("AAAAAA", "AAATAA"), choices)
#' matchBarcodes(c("AAAAAA", "AAATAA"), choices, substitutions=1)
#' matchBarcodes(c("AAAAAA", "AAATAA"), choices, reverse=TRUE)
#'
#' # Works with IUPAC codes in the barcodes:
#' choices <- c("AAARAA", "CCCYCC", "GGGMGG", "TTTSTT")
#' matchBarcodes(c("AAAAAA", "AAAGAA"), choices)
#'
#' @export
matchBarcodes <- function(sequences, choices, substitutions=0, reverse=FALSE) {
    out <- match_barcodes(sequences, choices, substitutions, reverse)
    DataFrame(index = out[[1]], mismatches = out[[2]])
}
