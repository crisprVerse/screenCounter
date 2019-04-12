#' Count single barcodes
#' 
#' Count the frequency of barcodes in a FASTQ file containing data for a single-end sequencing screen.
#'
#' @param fastq String containing the path to a FASTQ file containing single-end data,
#' or a connection object to such a file.
#' @param choices A \linkS4class{DataFrame} of sequences for the variable regions.
#' Each row should correspond to a barcode and each column should contain a character vector of sequences.
#' @param template A template for the barcode structure, see \code{\link{?createBarcodes}} for details.
#' @param substitutions Logical scalar specifying whether substitutions should be allowed when matching to variable regions.
#' @param deletions Logical scalar specifying whether deletions should be allowed when matching to variable regions.
#'
#' @details
#' In the simplest case, \code{choices} can be specified with a single column of variable sequences.
#' If \code{template=NULL}, the function will then directly use each of those sequences as the barcode with no constant regions.
#' More complex barcode structures are accommodated via \code{template}, e.g., with multiple variable regions and intervening and/or flanking constant spacers.
#' See \code{\link{createBarcodes}} for more details.
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
countSingleBarcodes <- function(fastq, choices, template=NULL, substitutions=TRUE, deletions=TRUE) {
    if (is.null(template)) {
        if (ncol(choices) > 1L) {
            stop("'template=NULL' only works with a single column of 'choices'")
        }
        template <- strrep("N", nchar(choices[,1]))
    }

    parsed <- parseBarcodeTemplate(template)
    n.pos <- parsed$variable$pos
    n.len <- parsed$variable$len
    constants <- parsed$constant

    # Validating 'choices'.
    nvariables <- length(n.pos)
    if (nvariables!=length(choices)) {
        stop("'length(choices)' is not equal to the number of stretches of N's")
    }
    for (i in seq_len(nvariables)) {
        if (!all(nchar(choices[[i]])==n.len[i])) {
            stop("each column of 'choices' must have same width as variable region in 'template'")
        }
    }

    # Choosing the C++ functions to use.
    if (nvariables==1L) {
        setupfun <- setup_barcodes_combo_solo
        countfun <- count_barcodes_combo_solo
        reportfun <- report_barcodes_combo_solo
    } else if (nvariables==2L) {
        setupfun <- setup_barcodes_combo_dual
        countfun <- count_barcodes_combo_dual
        reportfun <- report_barcodes_combo_dual
    } else {
        stop(sprintf("'ncol(choices)=%i' is not currently supported", nvariables))
    }

    ptr <- setup_barcodes_single(constants, as.list(barcodes)), substitutions, deletions)

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
