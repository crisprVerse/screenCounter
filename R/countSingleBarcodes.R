#' Count single barcodes
#' 
#' Count the frequency of barcodes in a FASTQ file containing data for a single-end sequencing screen.
#'
#' @param fastq String containing the path to a FASTQ file containing single-end data,
#' or a connection object to such a file.
#' @param choices A character vector of sequences for the variable regions, one per barcode.
#' @param flank5 String containing the constant sequence on the 5' flank of the variable region.
#' @param flank3 String containing the constant sequence on the 3' flank of the variable region.
#' @param template String containing the template for the barcode structure.
#' @param substitutions Integer scalar specifying the maximum number of substitutions when considering a match. 
#' @param insertions Integer scalar specifying the maximum number of insertions when considering a match. 
#' @param deletions Integer scalar specifying the maximum number of deletions when considering a match. 
#' @param total.edits Integer scalar specifying the maximum number of edits when considering a match.
#' @param strand String specifying which strand of the read to search.
#' @param files A character vector of paths to FASTQ files.
#' @param ... Further arguments to pass to \code{countSingleBarcodes}.
#' @param withDimnames A logical scalar indicating whether the rows and columns should be named.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how parallelization is to be performed across files.
#'
#' @details
#' If \code{template} is specified, it will be used to define the flanking regions.
#' Any user-supplied values of \code{flank5} and \code{flank3} will be ignored.
#' Note that, for this function, the template should only contain a single variable region.
#' See \code{\link{parseBarcodeTemplate}} for more details.
#'
#' We can handle sequencing errors across the entire barcode sequence (including variable and flanking regions) 
#' by setting \code{substitutions}, \code{deletions} and \code{insertions} to accept imperfect matches. 
#' If \code{total.edits} is specified, the total number of edits is capped regardless of the individual values for each edit type.
#' The default of \code{total.edits=2} means that only 2 edits are allowed for a match, even if \code{substitutions + insertions + deletions} is greater than 2.
#'
#' If \code{strand="both"}, the original read sequence will be searched first.
#' If no match is found, the sequence is reverse-complemented and searched again.
#' Other settings of \code{strand} will only search one or the other sequence.
#' The most appropriate choice depends on both the sequencing protocol and the design (i.e., position and length) of the barcode.
#'
#' @return 
#' \code{countSingleBarcodes} will return a \linkS4class{DataFrame} containing:
#' \itemize{
#' \item \code{choices}, a character vector equal to the input \code{choices}.
#' \item \code{counts}, an integer vector of length equal to \code{nrow(choices)} containing the frequency of each barcode.
#' }
#' The metadata contains \code{nreads}, an integer scalar containing the total number of reads in \code{fastq}.
#' 
#' \code{matrixOfSingleBarcodes} will return a \linkS4class{SummarizedExperiment} object containing:
#' \itemize{
#' \item An integer matrix named \code{"counts"}, where each column is the output of \code{countSingleBarcodes} for each file in \code{files}.
#' \item Row metadata containing a character vector \code{choices}, the sequence of the variable region of each barcode for each row.
#' \item Column metadata containing a character vector \code{files}, the path to each file;
#' an integer vector \code{nreads}, containing the total number of reads in each file;
#' and \code{nmapped}, containing the number of reads assigned to a barcode in the output count matrix.
#' }
#' If \code{withDimnames=TRUE}, row names are set to \code{choices} while column names are \code{basename(files)}.
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
#' matrixOfSingleBarcodes(c(tmp, tmp), choices=known.pool,
#'     flank5="CAGCTACGTACG", flank3="CCAGCTCGATCG")
#' @export
#' @importFrom ShortRead FastqStreamer yield sread
#' @importFrom S4Vectors DataFrame metadata<-
countSingleBarcodes <- function(fastq, choices, flank5, flank3, template=NULL, 
    substitutions=0, insertions=0, deletions=0, total.edits=2,
    strand=c("both", "original", "reverse")) 
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

    strand <- match.arg(strand)
    use.forward <- strand %in% c("original", "both")
    use.reverse <- strand %in% c("reverse", "both")

    ptr <- setup_barcodes_single(constants, choices, substitutions, insertions, deletions, total.edits)
    incoming <- FastqStreamer(fastq) 
    on.exit(close(incoming))

    N <- 0L
    all.available <- integer(length(choices))
    while (length(fq <- yield(incoming))) {
        current <- count_barcodes_single(sread(fq), ptr, use.forward, use.reverse)
        all.available <- all.available + current
        N <- N + length(fq)
    }

    out <- DataFrame(choices=choices, counts=all.available)
    metadata(out)$nreads <- N
    out
}

#' @rdname countSingleBarcodes
#' @export
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame metadata
matrixOfSingleBarcodes <- function(files, choices, ..., withDimnames=TRUE, BPPARAM=SerialParam()) {
    out <- bplapply(files, FUN=countSingleBarcodes, choices=choices, ..., BPPARAM=BPPARAM)
    mat <- do.call(cbind, lapply(out, "[[", "counts"))
    nreads <- vapply(out, function(x) metadata(x)$nreads, FUN.VALUE=0L)

    se <- SummarizedExperiment(list(counts=mat),
        rowData=DataFrame(choices=choices),
        colData=DataFrame(paths=files, nreads=nreads, nmapped=colSums(mat)))

    if (withDimnames) {
        rownames(se) <- choices
        colnames(se) <- basename(files)
    }
    se 
}
