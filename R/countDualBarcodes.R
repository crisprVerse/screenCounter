#' Count dual barcodes
#' 
#' Count the frequency of dual barcodes in a dataset for a paired-end sequencing screen.
#'
#' @param fastq Character vector of length 2, containing paths to two FASTQ files with paired-end data.
#' @param choices A \linkS4class{DataFrame} with two character columns specifying valid combinations of variable regions.
#' The first column contains sequences for barcode 1 while the second column contains sequences for barcode 2.
#' @param flank5 Character vector of length 2 containing the constant sequence on the 5' flank of the variable region for barcodes 1 and 2, respectively.
#' Alternatively, a string can be supplied if the constant sequence is the same for each barcode.
#' @param flank3 Character vector of length 2 containing the constant sequence on the 3' flank of the variable region for barcodes 1 and 2, respectively.
#' Alternatively, a string can be supplied if the constant sequence is the same for each barcode.
#' @param template Character vector of length 2 containing the template for the structure of barcodes 1 and 2, respectively.
#' Alternatively, a string can be supplied if the template is the same for each barcode.
#' @param substitutions Integer vector of length 2 specifying how many substitutions should be allowed for barcodes 1 and 2, respectively.
#' Alternatively, an integer scalar can be supplied if this is the same for each barcode.
#' @param find.best Logical scalar indicating whether to search each read for the best match.
#' Defaults to stopping at the first match.
#' @param strand Character vector of length 2 specifying which strand of the read to search (\code{"original"}, \code{"reverse"}) for each barcode.
#' Alternatively, a string can be supplied if this is the same for each barcode.
#' @param randomized Logical scalar indicating whether the first FASTQ file always contains the first barcode in \code{choices}.
#' If not, the opposite orientation is also searched.
#' @param include.invalid Logical scalar indicating whether counts for invalid barcode combinations should also be returned.
#' @param num.threads Integer scalar specifying the number of threads to use to process a single file.
#' @param files A list of character vectors of length 2 containing paths to paired FASTQ files.
#' @param ... Further arguments to pass to \code{countDualBarcodes}.
#' @param withDimnames A logical scalar indicating whether the rows and columns should be named.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how parallelization is to be performed across files.
#'
#' @details
#' In a dual barcode experiment, each read of a paired-end sequencing experiment contains one barcode.
#' The goal is to count the frequency of each combination of barcodes across the read pairs.
#' This differs from \code{\link{countComboBarcodes}} in that (i) only a subset of combinations are valid
#' and (ii) the two barcodes occur on different reads.
#'
#' The interpretation of the arguments for matching each barcode to reads is similar to that of \code{\link{countSingleBarcodes}}.
#' Each barcode in the combination can be associated with different search parameters;
#' for example, the search for the \dQuote{first} barcode in \code{choices[,1]}
#' will be performed with \code{flank5[1]}, \code{flank3[1]}, \code{substitutions[1]}, \code{strand[1]}, etc.
#' 
#' By default, the first FASTQ file is assumed to contain the first barcode (i.e., \code{choices[,1]})
#' while the second file is assumed to contain the second barcode (\code{choices[,2]}).
#' However, if \code{randomized=TRUE}, the orientation is assumed to be random such that
#' the first FASTQ file may contain the second barcode and so on.
#' In such cases, both orientations will be searched to identify a valid combination.
#'
#' We can handle sequencing errors by setting \code{substitutions} to a value greater than zero.
#' This will consider substitutions in both the variable region as well as the constant flanking regions.
#'
#' By default, the function will stop at the first match that satisfies the requirements above.
#' If \code{find.best=TRUE}, we will instead try to find the best match with the fewest mismatches.
#' If there are multiple matches with the same number of mismatches, the read is discarded to avoid problems with ambiguity.
#'
#' @return 
#' By default, \code{countDualBarcodes} will return \code{choices} with an additional \code{counts} column.
#' This is an integer vector of length equal to \code{nrow(choices)} containing the frequency of each barcode combination.
#' The metadata contains \code{npairs}, the total number of read pairs processed by the function.
#' 
#' \code{matrixOfDualBarcodes} will return a \linkS4class{SummarizedExperiment} object containing:
#' \itemize{
#' \item An integer matrix named \code{"counts"}, where each column is the output of \code{countDualBarcodes} for each file in \code{files}.
#' \item Row metadata containing a character vector \code{choices}, the sequences of the variable region of the two barcodes for each row.
#' \item Column metadata containing the character vectors \code{paths1} and \code{paths2}, storign the path to each pair of FASTQ files;
#' integer vectors corresponding to the metadata described above for \code{countDualBarcodes};
#' and \code{nmapped}, containing the number of read pairs assigned to a barcode combination in the output count matrix.
#' }
#' If \code{withDimnames=TRUE}, row names are set to \code{choices} while column names are \code{basename(files)}.
#'
#' If \code{include.invalid=TRUE}, each row contains all observed combinations in addition to those in \code{choices}.
#' The DataFrame (or \code{\link{rowData}} of the SummarizedExperiment) gains a \code{valid} field specifying if a combination is valid, i.e., present in \code{choices},
#' The metadata also gains the following fields:
#' \itemize{
#' \item \code{invalid.pair}, the number of read pairs with matches for each barcode but do not form a valid combination.
#' \item \code{barcode1.only}, the number of read pairs that only match to barcode 1.
#' \item \code{barcode2.only}, the number of read pairs that only match to barcode 2.
#' }
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
#' read1 <- sprintf("CAGCTACGTACG%sCCAGCTCGATCG",
#'    sample(known.pool1, N, replace=TRUE))
#' names(read1) <- seq_len(N)
#' 
#' read2 <- sprintf("TGGGCAGCGACA%sACACGAGGGTAT",
#'    sample(known.pool2, N, replace=TRUE))
#' names(read2) <- seq_len(N)
#' 
#' library(Biostrings)
#' tmp <- tempfile()
#' tmp1 <- paste0(tmp, "_1.fastq")
#' writeXStringSet(DNAStringSet(read1), filepath=tmp1, format="fastq")
#' tmp2 <- paste0(tmp, "_2.fastq")
#' writeXStringSet(DNAStringSet(read2), filepath=tmp2, format="fastq")
#'
#' # Counting the combinations.
#' countDualBarcodes(c(tmp1, tmp2), choices=choices, 
#'     template=c("CAGCTACGTACGNNNNNNNNNCCAGCTCGATCG",
#'                "TGGGCAGCGACANNNNNNNNNACACGAGGGTAT"))
#'
#' countDualBarcodes(c(tmp1, tmp2), choices=choices,
#'     flank5=c("CAGCTACGTACG", "TGGGCAGCGACA"),
#'     flank3=c("CCAGCTCGATCG", "ACACGAGGGTAT"))
#'
#' matrixOfDualBarcodes(list(c(tmp1, tmp2), c(tmp1, tmp2)),
#'     choices=choices,
#'     flank5=c("CAGCTACGTACG", "TGGGCAGCGACA"),
#'     flank3=c("CCAGCTCGATCG", "ACACGAGGGTAT"))
#' @export
#' @importFrom S4Vectors DataFrame metadata<- countMatches selfmatch
countDualBarcodes <- function(
    fastq, 
    choices, 
    flank5, 
    flank3, 
    template=NULL, 
    substitutions=0, 
    find.best=FALSE,
    strand="original", 
    randomized=FALSE, 
    include.invalid=FALSE, 
    num.threads=1)
{
    temp.out <- .create_paired_templates(template, flank5, flank3, choices)
    template1 <- temp.out[[1]]
    template2 <- temp.out[[2]]

    substitutions <- rep(substitutions, length.out=2)

    strand <- .verify_strand(strand)
    strand1 <- strand[1] == "reverse"
    strand2 <- strand[2] == "reverse"

    output <- count_dual_barcodes(
        fastq[1], template1, strand1, substitutions[1], as.character(choices[,1]),
        fastq[2], template2, strand2, substitutions[2], as.character(choices[,2]),
        randomized, !find.best, include.invalid, num.threads
    )

    if (!include.invalid) {
        choices$counts <- output[[1]]
        metadata(choices) <- list(npairs = output[[2]])
        return(choices)

    } else {
        combined <- .attach_invalid_counts(choices, valid.counts=output[[1]], invalid.combos=output[[2]][[1]], invalid.counts=output[[2]][[2]])
        metadata(combined) <- list(npairs = output[[3]], barcode1.only = output[[4]], barcode2.only = output[[5]], invalid.pair = sum(others$counts))
        return(combined)
    }
}

.create_paired_templates <- function(template, flank5, flank3, choices) {
    if (!is.null(template)) {
        template <- rep(template, length.out=2)
        template1 <- gsub("[nN]", "-", template[1])
        template2 <- gsub("[nN]", "-", template[2])
    } else {
        flank5 <- rep(flank5, length.out=2)
        flank3 <- rep(flank3, length.out=2)
        template1 <- paste0(flank5[1], strrep("-", nchar(as.character(choices[1,1]))), flank3[1])
        template2 <- paste0(flank5[2], strrep("-", nchar(as.character(choices[1,2]))), flank3[2])
    }
    c(template1, template2)
}

.verify_strand <- function(strand) {
    strand <- rep(strand, length.out=2)
    for (s in seq_along(strand)) {
        strand[s] <- match.arg(strand[s], c("original", "reverse"))
    }
    strand
}

.attach_invalid_counts <- function(choices, valid.counts, invalid.combos, invalid.counts) {
    choices$counts <- valid.counts
    choices$valid <- !logical(nrow(choices))

    others <- DataFrame(
        X1 = choices[invalid.combos[1,] + 1L,1], 
        X2 = choices[invalid.combos[2,] + 1L,2]
    )
    colnames(others) <- colnames(choices)[1:2]
    others$counts <- invalid.counts
    others$valid <- logical(nrow(others))

    if (!is.null(rownames(choices))) {
        rownames(others) <- character(nrow(others))
    }
    rbind(choices, others)
}

#' @export
#' @rdname countDualBarcodes
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom BiocParallel bplapply SerialParam
matrixOfDualBarcodes <- function(files, choices, ..., withDimnames=TRUE, include.invalid=FALSE, BPPARAM=SerialParam()) {
    out <- bplapply(files, FUN=countDualBarcodes, choices=choices, ..., include.invalid=include.invalid, BPPARAM=BPPARAM)
    se <- .post_process_dual_barcode_matrix(out, choices, include.invalid)

    colData(se) <- cbind(
        DataFrame(
            paths1=vapply(files, "[", i=1, ""), 
            paths2=vapply(files, "[", i=2, "")
        ),
        colData(se)
    )

    if (withDimnames) {
        colnames(se) <- basename(se$paths1)
    }
    se 

}

#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame metadata %in% metadata
.post_process_dual_barcode_matrix <- function(out, choices, include.invalid) {
    if (include.invalid) {
        tmp <- out
        for (i in seq_along(out)) {
            current <- out[[i]]
            df <- DataFrame(counts=current$counts)
            df$combinations <- current[,seq_len(ncol(current) - 1L)]
            tmp[[i]] <- df
        }

        combined <- do.call(combineComboCounts, tmp)
        mat <- combined$counts
        choices2 <- combined$combinations
        choices2$valid <- choices2 %in% choices
        metadata(choices2) <- list()
        choices <- choices2

    } else {
        mat <- do.call(cbind, lapply(out, "[[", "counts"))
    }

    fields <- names(metadata(out[[1]]))
    output <- list()
    for (f in fields) {
        output[[f]] <- vapply(out, function(x) metadata(x)[[f]], FUN.VALUE=0L)
    }

    SummarizedExperiment(list(counts=mat), rowData=choices, colData=output)
}
