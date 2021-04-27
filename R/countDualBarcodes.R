#' Count dual barcodes
#' 
#' Count the frequency of dual barcodes in a dataset for a paired-end sequencing screen.
#'
#' @param fastq Character vector of length two containing the path to two FASTQ files containing paired-end data.
#' Alternatively, a list of length two containing connection objects to such files.
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
#' @param insertions Integer vector of length 2 specifying how many insertions should be allowed for barcodes 1 and 2, respectively.
#' Alternatively, an integer scalar can be supplied if this is the same for each barcode.
#' @param deletions Integer vector of length 2 specifying how many deletions should be allowed for barcodes 1 and 2, respectively.
#' Alternatively, an integer scalar can be supplied if this is the same for each barcode.
#' @param total.edits Integer vector of length 2 specifying how many total edits should be allowed for barcodes 1 and 2, respectively.
#' Alternatively, an integer scalar can be supplied if this is the same for each barcode.
#' @param strand Character vector of length 2 specifying which strand of the read to search 
#' (\code{"both"}, \code{"original"}, \code{"reverse"}) for each barcode.
#' Alternatively, a string can be supplied if this is the same for each barcode.
#' @param randomized Logical scalar indicating whether the first FASTQ file always contains the first barcode in \code{choices}.
#' If not, the opposite orientation is also searched.
#' @param include.invalid Logical scalar indicating whether counts for invalid barcode combinations should also be returned.
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
#' The behavior of this function with respect to the actual matching of read sequences to each individual barcode 
#' is the same as that of \code{\link{countSingleBarcodes}}.
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
#' We can handle sequencing errors across the entire barcode sequence (including variable and flanking regions) 
#' by setting \code{substitutions}, \code{deletions} and \code{insertions} to accept imperfect matches. 
#' If \code{total.edits} is specified, the total number of edits is capped regardless of the individual values for each edit type.
#' The default of \code{total.edits=2} means that only 2 edits are allowed for a match, even if \code{substitutions + insertions + deletions} is greater than 2.
#'
#' @return 
#' By default, \code{countDualBarcodes} will return \code{choices} with an additional \code{counts} column
#' containing an integer vector of length equal to \code{nrow(choices)} containing the frequency of each barcode combination.
#' The metadata contains:
#' \itemize{
#' \item \code{none}, the number of read pairs with no matches to either barcode.
#' \item \code{barcode1.only}, the number of read pairs that only match to barcode 1.
#' \item \code{barcode2.only}, the number of read pairs that only match to barcode 2.
#' \item \code{invalid.pair}, the number of read pairs with matches for each barcode 
#' but do not form a valid combination in \code{choices}.
#' \item \code{provided.orientation}, the number of read pairs that match a barcode combination in the provided orientation,
#' i.e., the first and second FASTQ files contain the first and second barcodes, respectively.
#' Only reported when \code{randomized=TRUE}.
#' \item \code{other.orientation}, the number of read pairs that match a barcode combination in the other orientation,
#' i.e., the first and second FASTQ files contain the second and first barcodes, respectively.
#' Only reported when \code{randomized=TRUE}.
#' }
#' 
#' \code{matrixOfDualBarcodes} will return a \linkS4class{SummarizedExperiment} object containing:
#' \itemize{
#' \item An integer matrix named \code{"counts"}, where each column is the output of \code{countDualBarcodes} for each file in \code{files}.
#' \item Row metadata containing a character vector \code{choices}, the sequences of the variable region of the two barcodes for each row.
#' \item Column metadata containing a character vector \code{files}, the path to each file;
#' integer vectors corresponding to the metadata described above for \code{countDualBarcodes};
#' and \code{nmapped}, containing the number of read pairs assigned to a barcode combination in the output count matrix.
#' }
#' If \code{withDimnames=TRUE}, row names are set to \code{choices} while column names are \code{basename(files)}.
#'
#' If \code{include.invalid=TRUE}, each row contains all observed combinations in addition to those in \code{choices}.
#' The \code{\link{rowData}} contains an additional \code{valid} field specifying if a combination is valid, i.e., present in \code{choices}.
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
#' @importFrom BiocGenerics anyDuplicated match %in%
countDualBarcodes <- function(fastq, choices, flank5="", flank3="", template=NULL, 
    substitutions=0, insertions=0, deletions=0, total.edits=2,
    strand="original", randomized=FALSE, include.invalid=FALSE)
{
    # Checking the choices.
    condensed <- as.list(choices)
    condensed <- lapply(choices, unique)
    for (i in seq_along(choices)) {
        choices[[i]] <- match(choices[[i]], condensed[[i]])
    }
    if (anyDuplicated(choices)) {
        stop("'choices' must contain only unique combinations")
    }

    # Checking the search parameters.
    if (!is.null(template)) {
        template <- rep(template, length.out=2)
        flank5 <- flank3 <- character(0)

        for (i in seq_along(template)) {
            parsed <- parseBarcodeTemplate(template[i])
            constants <- parsed$constant
            if (length(constants)!=2L) {
                stop("template must contain exactly one variable region")
            }
            flank5 <- c(flank5, constants[1])
            flank3 <- c(flank3, constants[2])
        }
    } else {
        flank5 <- rep(flank5, length.out=2)
        flank3 <- rep(flank3, length.out=2)
    }

    substitutions <- rep(substitutions, length.out=2)
    insertions <- rep(insertions, length.out=2)
    deletions <- rep(deletions, length.out=2)
    total.edits <- rep(total.edits, length.out=2)
    strand <- rep(strand, length.out=2)

    # Performing the search.
    collected <- .dual_identifier(fastq, condensed, flank5=flank5, flank3=flank3, 
         substitutions=substitutions, insertions=insertions, deletions=deletions, total.edits=total.edits, 
         strand=strand, randomized=randomized)

    observed <- DataFrame(collected[[1]], collected[[2]])
    colnames(observed) <- colnames(choices)
    status <- (observed[,1] != 0L) + (observed[,2] != 0L) * 2L
    observed <- observed[status==3L,,drop=FALSE]

    # Generating counts.
    assignments <- selfmatch(observed)
    keep <- !duplicated(assignments)
    uniq.out <- observed[keep,,drop=FALSE]
    uniq.counts <- countMatches(assignments[keep], assignments)

    m <- match(choices, uniq.out)
    valid.counts <- uniq.counts[m]
    valid.counts[is.na(m)] <- 0L
    output <- choices
    output$counts <- valid.counts

    is.invalid <- setdiff(seq_len(nrow(uniq.out)), m)
    invalid.counts <- uniq.counts[is.invalid]

    # Appending on invalid pairs.
    if (include.invalid) {
        invalids <- uniq.out[is.invalid,,drop=FALSE]
        invalids$counts <- invalid.counts
        invalids$valid <- logical(nrow(invalids)) 
        if (!is.null(rownames(output))) {
            # Preserve the row names by giving something to combine with.
            rownames(invalids) <- rep("invalid", nrow(invalids))
        }

        output$valid <- !logical(nrow(output))
        output <- rbind(output, invalids)
    }

    # Restoring sequence identities.
    output[,1] <- condensed[[1]][output[,1]]
    output[,2] <- condensed[[2]][output[,2]]

    metadata(output) <- c(
        list(
            none=sum(status==0L), 
            barcode1.only=sum(status==1L),
            barcode2.only=sum(status==2L), 
            invalid.pair=sum(invalid.counts),
            nmapped=sum(valid.counts)
        ),
        collected$diagnostics
    )

    output
}

#' @importFrom ShortRead FastqStreamer yield sread
.dual_identifier <- function(fastq, condensed, flank5, flank3, substitutions, insertions, deletions, total.edits, strand, randomized) {
    ptr1 <- setup_barcodes_dual(c(flank5[1], flank3[1]), condensed[[1]], substitutions[1], insertions[1], deletions[1], total.edits[1])
    ptr2 <- setup_barcodes_dual(c(flank5[2], flank3[2]), condensed[[2]], substitutions[2], insertions[2], deletions[2], total.edits[2])

    incoming1 <- FastqStreamer(fastq[[1]]) 
    on.exit(close(incoming1))
    incoming2 <- FastqStreamer(fastq[[2]]) 
    on.exit(close(incoming2), add=TRUE)

    for (i in seq_along(strand)) {
        strand[i] <- match.arg(strand[i], c("original", "reverse"))
    }

    collected1 <- collected2 <- list()
    diagnostics <- c(0L, 0L)
    counter <- 1L

    repeat {
        fq1 <- yield(incoming1)
        fq2 <- yield(incoming2)
        if (length(fq1)==0 && length(fq2)==0) {
            break
        }

        output <- count_barcodes_dual(sread(fq1), sread(fq2), ptr1, ptr2, forward1=(strand[1]=="original"), forward2=(strand[2]=="original"), randomized=randomized)
        collected1[[counter]] <- output[[1]]
        collected2[[counter]] <- output[[2]]
        diagnostics <- diagnostics + output[[3]]
        counter <- counter + 1L
    }

    output <- list(first=unlist(collected1), second=unlist(collected2))
    if (randomized) {
        output$diagnostics <- list(provided.orientation=diagnostics[1], other.orientation=diagnostics[2])
    }
    output
}

#' @rdname countDualBarcodes
#' @export
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame metadata
matrixOfDualBarcodes <- function(files, choices, ..., withDimnames=TRUE, include.invalid=FALSE, BPPARAM=SerialParam()) {
    out <- bplapply(files, FUN=countDualBarcodes, choices=choices, ..., include.invalid=include.invalid, BPPARAM=BPPARAM)

    if (include.invalid) {
        tmp <- out
        for (i in seq_along(out)) {
            current <- out[[i]]
            df <- DataFrame(counts=current$counts)
            df$combinations <- current[,1:2]
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

    se <- SummarizedExperiment(list(counts=mat), rowData=choices,
        colData=DataFrame(paths1=vapply(files, "[", i=1, ""), 
            paths2=vapply(files, "[", i=2, ""), output))

    if (withDimnames) {
        colnames(se) <- basename(se$paths1)
    }
    se 
}
