#' Count single barcodes
#' 
#' Count the frequency of barcodes in a FASTQ file containing data for a single-end sequencing screen.
#'
#' @param fastq String containing the path to a FASTQ file containing single-end data,
#' or a connection object to such a file.
#' @param choices A character vector of sequences for the variable regions, one per barcode.
#' @param flank5 String containing the constant sequence on the 5' flank of the variable region.
#' @param flank3 String containing the constant sequence on the 3' flank of the variable region.
#' @param template A template for the barcode structure, see \code{\link{?createBarcodes}} for details.
#' @param substitutions Logical scalar specifying whether substitutions should be allowed when matching to variable regions.
#' @param deletions Logical scalar specifying whether deletions should be allowed when matching to variable regions.
#'
#' @details
#' If \code{template} is specified, it will be used to define the flanking regions.
#' Any user-supplied values of \code{flank5} and \code{flank3} will be ignored.
#' Note that, for this function, the template should only contain a single variable region.
#' See \code{\link{parseBarcodeTemplate}} for more details.
#'
#' If \code{substitutions=TRUE}, only one mismatch is allowed across all variable regions,
#' \emph{not} per variable region.
#' Similarly, if \code{deletions=TRUE}, only one deletion is allowed across all variable regions.
#' If both are set, only one deletion or mismatch is allowed across all variable regions,
#' i.e., there is a maximum edit distance of 1 from any possible reference combination.
#'
#' @return An integer vector of length equal to \code{nrow(choices)},
#' containing the frequency of each barcode.
#' This is named with the row names of \code{choices}.
#'
#' @author Aaron Lun
#' @examples
#' # Creating an example dual barcode sequencing experiment.
#' known.pool <- c("AGAGAGAGA", "CTCTCTCTC",
#'     "GTGTGTGTG", "CACACACAC")
#' 
#' N <- 1000
#' barcodes <- sprintf("CAGCTACGTACG%sCCAGCTCGATCG",
#'    sample(known.pool, N, replace=TRUE))
#' names(barcodes) <- seq_len(N)
#' 
#' library(Biostrings)
#' tmp <- tempfile(fileext=".fastq")
#' writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")
#'
#' # Counting the combinations.
#' countSingleBarcodes(tmp, choices=known.pool,
#'     template="CAGCTACGTACGNNNNNNNNNCCAGCTCGATCG")
#'
#' countSingleBarcodes(tmp, choices=known.pool,
#'     flank5="CAGCTACGTACG", flank3="CCAGCTCGATCG")
#'
#' @export
#' @importFrom ShortRead FastqStreamer yield sread
countSingleBarcodes <- function(fastq, choices, flank5, flank3, 
    template=NULL, substitutions=TRUE, deletions=TRUE) 
{
    if (!is.null(template)) {
        parsed <- parseBarcodeTemplate(template)
        constants <- parsed$constant
        if (length(constants)!=2L) {
            stop("template must contain exactly one variable region")
        }
    } else {
        constants <- as.character(c(flank5, flank3))
    }

    ptr <- setup_barcodes_single(constants, list(choices), substitutions, deletions)
    incoming <- FastqStreamer(fastq) 
    on.exit(close(incoming))
    while (length(fq <- yield(incoming))) {
        seqs <- as.character(sread(fq))
        count_barcodes_single(seqs, ptr)
    }

    all.available <- report_barcodes_single(ptr)
    names(all.available) <- rownames(choices)
    all.available
}
