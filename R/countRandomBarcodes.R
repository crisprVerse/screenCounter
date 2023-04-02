#' Count random barcodes
#' 
#' Count the frequency of random barcodes in a FASTQ file containing data for a single-end sequencing screen.
#' This differs from \code{\link{countSingleBarcodes}} in that the barcode is completely random rather than being drawn from a known pool of sequences.
#'
#' @param template String containing the template for the barcode structure.
#' See \code{\link{parseBarcodeTemplate}} for more details.
#' @inheritParams countSingleBarcodes
#'
#' @details
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
#' \code{countRandomBarcodes} will return a \linkS4class{DataFrame} containing:
#' \itemize{
#' \item \code{sequences}, a character vector containing the sequences of the random barcodes in the variable region.
#' \item \code{counts}, an integer vector containing the frequency of each barcode.
#' }
#' The metadata contains \code{nreads}, an integer scalar containing the total number of reads in \code{fastq}.
#' 
#' \code{matrixOfRandomBarcodes} will return a \linkS4class{SummarizedExperiment} object containing:
#' \itemize{
#' \item An integer matrix named \code{"counts"}, where each column is the output of \code{countRandomBarcodes} for each file in \code{files}.
#' \item Row metadata containing a character vector \code{sequences}, the sequence of the variable region of each barcode for each row.
#' \item Column metadata containing a character vector \code{files}, the path to each file;
#' an integer vector \code{nreads}, containing the total number of reads in each file;
#' and \code{nmapped}, containing the number of reads assigned to a barcode in the output count matrix.
#' }
#' If \code{withDimnames=TRUE}, row names are set to \code{sequences} while column names are \code{basename(files)}.
#'
#' @author Aaron Lun
#' @examples
#' # Creating an example dataset.
#' N <- 1000
#' randomized <- lapply(1:N, function(i) {
#'     paste(sample(c("A", "C", "G", "T"), 8, replace=TRUE), collapse="")
#' })
#' barcodes <- sprintf("CAGCTACGTACG%sCCAGCTCGATCG", randomized)
#' names(barcodes) <- seq_len(N)
#' 
#' library(Biostrings)
#' tmp <- tempfile(fileext=".fastq")
#' writeXStringSet(DNAStringSet(barcodes), filepath=tmp, format="fastq")
#'
#' # Counting the sequences:
#' countRandomBarcodes(tmp, template="CAGCTACGTACGNNNNNNNNCCAGCTCGATCG")
#'
#' matrixOfRandomBarcodes(c(tmp, tmp), template="CAGCTACGTACGNNNNNNNNCCAGCTCGATCG")
#' @export
#' @importFrom S4Vectors DataFrame metadata<-
countRandomBarcodes <- function(
    fastq, 
    template,
    substitutions=0, 
    find.best=FALSE,
    strand=c("both", "original", "reverse"),
    num.threads = 1) 
{
    template <- gsub("N", "-", template)
    strand <- c(original = 0L, reverse = 1L, both = 2L)[match.arg(strand)]
    results <- count_random_barcodes(fastq, template, strand, substitutions, !find.best, num.threads)

    o <- order(results[[1]][[1]])
    out <- DataFrame(sequences=results[[1]][[1]][o], counts=results[[1]][[2]][o])
    metadata(out)$nreads <- results[[2]]
    out
}

#' @rdname countRandomBarcodes
#' @export
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame metadata
matrixOfRandomBarcodes <- function(files, ..., withDimnames=TRUE, BPPARAM=SerialParam()) {
    out <- bplapply(files, FUN=countRandomBarcodes, ..., BPPARAM=BPPARAM)

    all.seq <- sort(Reduce(union, lapply(out, "[[", "sequences")))
    mat <- matrix(0L, length(all.seq), length(out))
    for (i in seq_along(out)) {
        m <- match(out[[i]]$sequences, all.seq)
        mat[m,i] <- out[[i]]$counts
    }

    nreads <- vapply(out, function(x) metadata(x)$nreads, FUN.VALUE=0L)

    se <- SummarizedExperiment(list(counts=mat),
        rowData=DataFrame(sequences=all.seq),
        colData=DataFrame(paths=files, nreads=nreads, nmapped=colSums(mat)))

    if (withDimnames) {
        rownames(se) <- all.seq
        colnames(se) <- basename(files)
    }
    se 
}
