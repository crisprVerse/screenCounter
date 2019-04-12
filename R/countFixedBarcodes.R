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
#' countFixedBarcodes(tmp, choices=DataFrame(known.pool))
#'
#' @export
#' @importFrom ShortRead FastqStreamer yield sread
countFixedBarcodes <- function(fastq, choices, template=NULL, substitutions=TRUE, deletions=TRUE) {
    if (is.null(template)) {
        if (ncol(choices) > 1L) {
            stop("'template=NULL' only works with a single column of 'choices'")
        }
        template <- strrep("N", nchar(choices[1,1]))
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
        setupfun <- setup_barcodes_fixed_solo
        countfun <- count_barcodes_fixed_solo
        reportfun <- report_barcodes_fixed_solo
    } else if (nvariables==2L) {
        setupfun <- setup_barcodes_fixed_dual
        countfun <- count_barcodes_fixed_dual
        reportfun <- report_barcodes_fixed_dual
    } else {
        stop(sprintf("'ncol(choices)=%i' is not currently supported", nvariables))
    }

    ptr <- setupfun(constants, as.list(choices), substitutions, deletions)

    incoming <- FastqStreamer(fastq) 
    on.exit(close(incoming))
    while (length(fq <- yield(incoming))) {
        seqs <- as.character(sread(fq))
        countfun(seqs, ptr)
    }

    all.available <- reportfun(ptr)
    names(all.available) <- rownames(choices)
    all.available
}
