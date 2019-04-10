#' Count combinatorial barcodes
#'
#' Count combinatorial barcodes for screen sequencing experiments where entities are distinguished based on random combinations of a small pool of known sequences. 
#' 
#' @param fastq String containing the path to a FASTQ file containing single-end data,
#' or a connection object to such a file.
#' @param template A template for the barcode structure, see \code{\link{?createBarcodes}} for details.
#' @param choices A \linkS4class{List} of potential sequences for each variable region in \code{template}.
#' Each row should correspond to a barcode and each column should contain a character vector of sequences.
#' 
#' @details
#' Certain screen sequencing experiments take advantage of combinatorial complexity to generate a very large pool of unique barcode sequences.
#' Only a subset of all possible combinatorial barcodes will be used in any given experiment.
#' This function only counts the combinations that are actually observed, improving efficiency over a more conventional approach (i.e., to generate all possible combinations and use \code{\link{countSingleBarcodes}} to count their frequency).
#' 
#' @return A \linkS4class{DataFrame} where each row corresponds to a combinatorial barcode.
#' It contains \code{keys}, a nested \linkS4class{DataFrame} where each column corresponds to an element of \code{choices} and contains the indices of the sequences in each combinatorial barcode;
#' and \code{counts}, an integer vector containing the frequency of each barcode.
#'
#' @author Aaron Lun
#' @examples
#' # Creating an example dual barcode sequencing experiment.
#' library(Biostrings)
#' known.pool <- c("AGAGAGAGA", "CTCTCTCTC",
#'     "GTGTGTGTG", "CACACACAC")
#' 
#' N <- 1000
#' barcodes <- sprintf("ACGT%sACGT%sACGT",
#'    sample(known.pool, N, replace=TRUE),
#'    sample(known.pool, N, replace=TRUE))
#' names(barcodes) <- seq_len(N)
#' 
#' tmp <- tempfile(fileext=".fastq")
#' writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")
#'
#' # Counting the combinations.
#' output <- countComboBarcodes(tmp,
#'     template="ACGTNNNNNNNNNACGTNNNNNNNNNACGT",
#'     choices=list(first=known.pool, second=known.pool))
#' output$combination
#' head(output$count)
#'
#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom ShortRead FastqStreamer
countComboBarcodes <- function(fastq, template, choices) {
    positions <- .split_template(template)
    n.pos <- positions$pos
    n.len <- positions$len

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

    # Extracting constant regions.
    last.pos <- 1L
    constants <- character(length(n.pos)+1L)
    for (i in seq_along(n.pos)) {
        constants[i] <- substring(template, last.pos, n.pos[i]-1)
        last.pos <- n.pos[i] + n.len[i]
    }
    constants[length(n.pos)+1] <- substring(template, last.pos, nchar(template))

    # Choosing the C++ functions to use.
    if (nvariables==2L) {
        setupfun <- setup_barcodes_combo_dual
        countfun <- count_barcodes_combo_dual
        reportfun <- report_barcodes_combo_dual
    } else {
        stop(sprintf("'ncol(choices)=%i' is not currently supported", nvariables))
    }

    # Counting all pairs of barcodes. 
    ptr <- setupfun(constants, as.list(choices))

    incoming <- FastqStreamer(fastq) 
    on.exit(close(incoming))
    while (length(fq <- yield(incoming))) {
        seqs <- as.character(sread(fq))
        countfun(seqs, ptr)
    }

    output <- reportfun(ptr)
    keys <- do.call(DataFrame, output[[2]])
    colnames(keys) <- names(choices)
    DataFrame(combination=I(keys), count=output[[1]])
}
