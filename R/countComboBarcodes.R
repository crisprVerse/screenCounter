#' Count combinatorial barcodes
#'
#' Count combinatorial barcodes for screen sequencing experiments where entities are distinguished based on random combinations of a small pool of known sequences. 
#' 
#' @param fastq String containing the path to a FASTQ file containing single-end data,
#' or a connection object to such a file.
#' @param template A template for the barcode structure, see \code{?\link{parseBarcodeTemplate}} for details.
#' @param choices A \linkS4class{List} of character vectors, one per variable region in \code{template}.
#' Each vector should contain the potential sequences for the corresponding variable region.
#' @param substitutions Logical scalar specifying whether substitutions should be allowed when matching to variable regions.
#' @param deletions Logical scalar specifying whether deletions should be allowed when matching to variable regions.
#' @param strand String specifying which strand of the read to search.
#' @param indices Logical scalar indicating whether integer indices should be used to define each combinational barcode.
#' @param files A character vector of paths to FASTQ files.
#' @param ... Further arguments to pass to \code{countComboBarcodes}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how parallelization is to be performed across files.
#' 
#' @details
#' Certain screen sequencing experiments take advantage of combinatorial complexity to generate a very large pool of unique barcode sequences.
#' Only a subset of all possible combinatorial barcodes will be used in any given experiment.
#' This function only counts the combinations that are actually observed, improving efficiency over a more conventional approach (i.e., to generate all possible combinations and use \code{\link{countSingleBarcodes}} to count their frequency).
#'
#' If \code{substitutions=TRUE}, only one mismatch is allowed across all variable regions,
#' \emph{not} per variable region.
#' Similarly, if \code{deletions=TRUE}, only one deletion is allowed across all variable regions.
#' If both are set, only one deletion or mismatch is allowed across all variable regions,
#' i.e., there is a maximum edit distance of 1 from any possible reference combination.
#' 
#' If \code{strand="both"}, the original read sequence will be searched first.
#' If no match is found, the sequence is reverse-complemented and searched again.
#' Other settings of \code{strand} will only search one or the other sequence.
#' The most appropriate choice depends on both the sequencing protocol and the design (i.e., position and length) of the barcode.
#'
#' @return 
#' \code{countComboBarcodes} returns a \linkS4class{DataFrame} where each row corresponds to a combinatorial barcode.
#' It contains \code{keys}, a nested \linkS4class{DataFrame} that contains the sequences that define each combinatorial barcode;
#' and \code{counts}, an integer vector containing the frequency of each barcode.
#'
#' \code{matrixOfComboBarcodes} returns the same output as \code{\link{combineComboCounts}},
#' i.e., a \linkS4class{DataFrame} containing \code{keys} and a count matrix of frequencies for each combinatorial barcode in each \code{files}.
#' Columns of the count matrix are named with \code{files}.
#'
#' Each column of \code{keys} corresponds to a single variable region in \code{template} and thus one entry of \code{choices}.
#' By default, the sequences are reported directly as character vectors.
#' If \code{indices=FALSE}, each column contains the indices of the sequences in the corresponding entry of \code{choices}.
#' 
#' @author Aaron Lun
#' @examples
#' # Creating an example dual barcode sequencing experiment.
#' known.pool <- c("AGAGAGAGA", "CTCTCTCTC",
#'     "GTGTGTGTG", "CACACACAC")
#' 
#' N <- 1000
#' barcodes <- sprintf("ACGT%sACGT%sACGT",
#'    sample(known.pool, N, replace=TRUE),
#'    sample(known.pool, N, replace=TRUE))
#' names(barcodes) <- seq_len(N)
#' 
#' library(Biostrings)
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
#' matrixOfComboBarcodes(c(tmp, tmp),
#'     template="ACGTNNNNNNNNNACGTNNNNNNNNNACGT",
#'     choices=list(first=known.pool, second=known.pool))
#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom ShortRead FastqStreamer sread yield
countComboBarcodes <- function(fastq, template, choices, substitutions=FALSE, deletions=FALSE,
    strand=c("original", "reverse", "both"), indices=FALSE)
{
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
    if (nvariables==2L) {
        setupfun <- setup_barcodes_combo_dual
        countfun <- count_barcodes_combo_dual
        reportfun <- report_barcodes_combo_dual
    } else {
        stop(sprintf("'ncol(choices)=%i' is not currently supported", nvariables))
    }

    strand <- match.arg(strand)
    use.forward <- strand %in% c("original", "both")
    use.reverse <- strand %in% c("reverse", "both")

    # Counting all pairs of barcodes. 
    ptr <- setupfun(constants, as.list(choices), substitutions, deletions)

    incoming <- FastqStreamer(fastq) 
    on.exit(close(incoming))
    while (length(fq <- yield(incoming))) {
        countfun(sread(fq), ptr, use.forward, use.reverse)
    }

    output <- reportfun(ptr)
    keys <- do.call(DataFrame, output[[2]])
    if (!is.null(names(choices))) {
        colnames(keys) <- names(choices)
    }
    if (!indices) {
        for (i in seq_len(ncol(keys))) {
            keys[,i] <- choices[[i]][keys[,i]]
        }
    }
    DataFrame(combinations=I(keys), counts=output[[1]])
}

#' @rdname countComboBarcodes
#' @export
#' @importFrom BiocParallel SerialParam bplapply
matrixOfComboBarcodes <- function(files, ..., BPPARAM=SerialParam()) {
    out <- bplapply(files, FUN=countComboBarcodes, ..., BPPARAM=BPPARAM)
    names(out) <- files
    do.call(combineComboCounts, out)
}
