#' Count single barcodes
#' 
#' Count the frequency of barcodes in a FASTQ file containing data for a single-end sequencing screen.
#'
#' @param fastq String containing the path to a FASTQ file containing single-end data,
#' or a connection object to such a file.
#' @param choices A \linkS4class{DataFrame} of potential barcodes.
#' Each row should correspond to a barcode and each column should contain a character vector of sequences.
#' @param template A template for the barcode structure, see \code{\link{?createBarcodes}} for details.
#'
#' @details
#' In the simplest case, only \code{choices} needs to be specified with the barcode sequences.
#' However, this function can accommodate more complex barcode structures via \code{template},
#' e.g., with multiple variable regions and intervening spacers.
#' See \code{\link{createBarcodes}} for more details.
#'
#' @return An integer vector of length equal to \code{nrow(choices)},
#' containing the frequency of each barcode.
#' This is named with the row names of \code{choices}.
#'
#' @author Aaron Lun
#' @examples
#' # Re-using ShortRead's examples:
#' library(ShortRead)
#' sp <- SolexaPath(system.file('extdata', package='ShortRead'))
#' fl <- file.path(analysisPath(sp), "s_1_sequence.txt")
#'
#' countSingleBarcodes(fl,
#'    choices=DataFrame(c("AAA", "GGG", "CCC", "TTT")))
#'
#' @export
#' @importFrom ShortRead FastqStreamer yield sread
countSingleBarcodes <- function(fastq, choices, template="N") {
    barcodes <- createBarcodes(template, choices)
    ptr <- setup_barcodes_single(barcodes)

    incoming <- FastqStreamer(fastq) 
    on.exit(close(incoming))
    all.available <- integer(length(barcodes))

    while (length(fq <- yield(incoming))) {
        seqs <- as.character(sread(fq))
        output <- count_barcodes_single(seqs, ptr)
        output <- output[output >= 0L] 
        all.available <- all.available + tabulate(output, length(all.available))
    }

    names(all.available) <- names(barcodes)
    all.available
}
