#' Count dual barcodes in single-end data
#' 
#' Count the frequency of dual barcodes in a single-end sequencing screen.
#'
#' @param fastq Character vector containing a path to a FASTQ file.
#' @param choices A \linkS4class{DataFrame} with one or more character columns specifying valid combinations of variable regions.
#' Each column contains sequences for successive barcodes in \code{template}.
#' @param template String containing the template for the barcode structure.
#' The number of variable regions should be equal to the number of columns of \code{choices}.
#' @param substitutions Integer specifying how many substitutions should be allowed. 
#' @param find.best Logical scalar indicating whether to search each read for the best match.
#' Defaults to stopping at the first match.
#' @param strand String specifying the strand of the read to search (\code{"original"}, \code{"reverse"}).
#' @param include.invalid Logical scalar indicating whether counts for invalid barcode combinations should also be returned.
#' This is currently only enabled for \code{template} with 2 variable regions.
#' @param num.threads Integer scalar specifying the number of threads to use to process a single file.
#' @param files Character vectors containing paths to FASTQ files.
#' @param ... Further arguments to pass to \code{countDualBarcodesSingleEnd}.
#' @param withDimnames A logical scalar indicating whether the rows and columns should be named.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how parallelization is to be performed across files.
#'
#' @details
#' In a dual barcode experiment, each read of a single-end sequencing experiment contains a barcode element with multiple variable regions.
#' The goal is to count the frequency of each combination of barcodes.
#' However, unlike \code{\link{countComboBarcodes}}, only a subset of combinations are valid here as defined in \code{choices}.
#'
#' The interpretation of the arguments for matching each barcode to reads is similar to that of \code{\link{countSingleBarcodes}}.
#' The strand of the read to search is defined with \code{strand}, defaulting to searching both strands.
#' We can handle sequencing errors by setting \code{substitutions} to a value greater than zero.
#' This will consider substitutions in both the variable region as well as the constant flanking regions.
#'
#' By default, the function will stop at the first match that satisfies the requirements above.
#' If \code{find.best=TRUE}, we will instead try to find the best match with the fewest mismatches.
#' If there are multiple matches with the same number of mismatches, the read is discarded to avoid problems with ambiguity.
#'
#' @return 
#' By default, \code{countDualBarcodesSingleEnd} will return \code{choices} with an additional \code{counts} column.
#' This is an integer vector of length equal to \code{nrow(choices)} containing the frequency of each barcode combination.
#' The metadata contains \code{nreads}, the total number of reads processed by the function.
#' 
#' \code{matrixOfDualBarcodesSingleEnd} will return a \linkS4class{SummarizedExperiment} object containing:
#' \itemize{
#' \item An integer matrix named \code{"counts"}, where each column is the output of \code{countDualBarcodes} for each file in \code{files}.
#' \item Row metadata containing a character vector \code{choices}, the sequences of the variable region of the two barcodes for each row.
#' \item Column metadata containing a character vector \code{paths}, the path to each FASTQ file;
#' and integer vectors corresponding to the metadata described above for \code{countDualBarcodesSingleEnd}.
#' }
#' If \code{withDimnames=TRUE}, row names are set to \code{choices} while column names are \code{basename(files)}.
#'
#' If \code{include.invalid=TRUE}, each row contains all observed combinations in addition to those in \code{choices}.
#' The DataFrame (or \code{\link{rowData}} of the SummarizedExperiment) gains a \code{valid} field specifying if a combination is valid, i.e., present in \code{choices}.
#' The metadata also gains the \code{invalid.reads} field, containing the number of reads with matches for each barcode but do not form a valid combination.
#'
#' @author Aaron Lun
#' @examples
#' # Creating an example dual barcode sequencing experiment.
#' known.pool1 <- c("AGAGAGAGA", "CTCTCTCTC",
#'     "GTGTGTGTG", "CACACACAC")
#' known.pool2 <- c("ATATATATA", "CGCGCGCGC",
#'     "GAGAGAGAG", "CTCTCTCTC")
#' choices <- expand.grid(known.pool1, known.pool2)
#' choices <- DataFrame(barcode1=choices[,1], barcode2=choices[,2])
#' 
#' N <- 1000
#' read <- sprintf(
#'    "CAGCTACGTACG%sCCAGCTCGATCG%sACACGAGGGTAT",
#'    sample(known.pool1, N, replace=TRUE),
#'    sample(known.pool2, N, replace=TRUE)
#' )
#' names(read) <- seq_len(N)
#' 
#' library(Biostrings)
#' tmp <- tempfile(fileext=".fastq")
#' writeXStringSet(DNAStringSet(read), filepath=tmp, format="fastq")
#'
#' # Counting the combinations.
#' countDualBarcodesSingleEnd(tmp, choices=choices, 
#'     template="CAGCTACGTACGNNNNNNNNNCCAGCTCGATCGNNNNNNNNNACACGAGGGTAT")
#'
#' matrixOfDualBarcodesSingleEnd(c(tmp, tmp),
#'     choices=choices,
#'     template="CAGCTACGTACGNNNNNNNNNCCAGCTCGATCGNNNNNNNNNACACGAGGGTAT")
#' @export
#' @importFrom S4Vectors DataFrame metadata<- countMatches selfmatch
countDualBarcodesSingleEnd <- function(
    fastq, 
    choices, 
    template,
    substitutions=0, 
    find.best=FALSE,
    strand=c("both", "original", "reverse"),
    include.invalid=FALSE, 
    num.threads=1)
{
    template <- gsub("[nN]", "-", template)
    strand <- c(original=0L, reverse=1L, both=2L)[match.arg(strand)]

    output <- count_dual_barcodes_single_end(
        fastq, 
        template, 
        strand=strand, 
        mismatches=substitutions, 
        pools=as.list(choices), 
        use_first=!find.best, 
        diagnostics=include.invalid, 
        nthreads=num.threads
    )

    if (!include.invalid) {
        choices$counts <- output[[1]]
        metadata(choices) <- list(nreads = output[[2]])
        return(choices)

    } else {
        others <- .create_invalid_df(choices, invalid.combos=output[[2]][[1]], invalid.counts=output[[2]][[2]])
        choices$counts <- output[[1]]
        choices$valid <- !logical(nrow(choices))
        combined <- rbind(choices, others)
        metadata(combined) <- list(nreads = output[[3]], invalid.reads = sum(others$counts))
        return(combined)
    }
}

#' @export
#' @rdname countDualBarcodesSingleEnd
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom BiocParallel bplapply SerialParam
matrixOfDualBarcodesSingleEnd <- function(files, choices, ..., withDimnames=TRUE, include.invalid=FALSE, BPPARAM=SerialParam()) {
    out <- bplapply(files, FUN=countDualBarcodesSingleEnd, choices=choices, ..., include.invalid=include.invalid, BPPARAM=BPPARAM)
    se <- .post_process_dual_barcode_matrix(out, choices, include.invalid)
    colData(se) <- cbind(DataFrame(paths=files), colData(se))
    if (withDimnames) {
        colnames(se) <- basename(se$paths)
    }
    se 
}
