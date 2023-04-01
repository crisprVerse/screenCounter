#' Count single barcodes
#' 
#' Count the frequency of barcodes in a FASTQ file containing data for a single-end sequencing screen.
#'
#' @param fastq String containing the path to a FASTQ file containing single-end data.
#' @param choices A character vector of sequences for the variable regions, one per barcode.
#' @param flank5 String containing the constant sequence on the 5' flank of the variable region.
#' @param flank3 String containing the constant sequence on the 3' flank of the variable region.
#' @param template String containing the template for the barcode structure.
#' @param substitutions Integer scalar specifying the maximum number of substitutions when considering a match. 
#' @param find.best Logical scalar indicating whether to search each read for the best match.
#' Defaults to stopping at the first match.
#' @param strand String specifying which strand of the read to search.
#' @param num.threads Integer scalar specifying the number of threads to use to process a single file.
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
#' If \code{strand="both"}, the original read sequence will be searched first.
#' If no match is found, the sequence is reverse-complemented and searched again.
#' Other settings of \code{strand} will only search one or the other sequence.
#' The most appropriate choice depends on both the sequencing protocol and the design (i.e., position and length) of the barcode.
#'
#' We can handle sequencing errors by setting \code{substitutions} to a value greater than zero.
#' This will consider substitutions in both the variable region as well as the constant flanking regions.
#'
#' By default, the function will stop at the first match that satisfies the requirements above.
#' If \code{find.best=TRUE}, we will instead try to find the best match with the fewest mismatches.
#' If there are multiple matches with the same number of mismatches, the read is discarded to avoid problems with ambiguity.
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
#' @importFrom S4Vectors DataFrame metadata<-
countSingleBarcodes <- function(
    fastq, 
    choices, 
    flank5, 
    flank3, 
    template=NULL, 
    substitutions=0, 
    find.best=FALSE,
    strand=c("both", "original", "reverse"),
    num.threads = 1) 
{
    if (!is.null(template)) {
        template <- gsub("N", "-", template)
    } else {
        template <- paste0(flank5, strrep("-", nchar(choices)[1]), flank3)
    }

    strand <- c(original = 0L, reverse = 1L, both = 2L)[match.arg(strand)]
    results <- count_single_barcodes(fastq, template, strand, choices, substitutions, !find.best, num.threads)

    out <- DataFrame(choices=choices, counts=results[[1]])
    metadata(out)$nreads <- results[[2]]
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
