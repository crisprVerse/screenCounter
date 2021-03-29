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
#' (If both orientations yield a different valid combination, no count will be assigned.)
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
#' \item \code{original.orientation}, the number of read pairs that match a barcode combination in the original strand orientation.
#' Only reported when \code{randomized=TRUE}.
#' \item \code{reverse.orientation}, the number of read pairs that match a barcode combination in the reverse strand orientation.
#' Only reported when \code{randomized=TRUE}.
#' \item \code{ambiguous.orientation}, the number of read pairs that match a different barcode combination in each strand orientation.
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
    original <- choices
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

    # Searching two sequences for the matches.
    r1 <- .solo_identifier(fastq[[1]], choices=condensed[[1]], flank5=flank5[1], flank3=flank3[1],
        substitutions=substitutions[1], insertions=insertions[1], deletions=deletions[1], total.edits=total.edits[1],
        strand=strand[1])
    r2 <- .solo_identifier(fastq[[2]], choices=condensed[[2]], flank5=flank5[2], flank3=flank3[2],
        substitutions=substitutions[2], insertions=insertions[2], deletions=deletions[2], total.edits=total.edits[2],
        strand=strand[2])

    observed <- DataFrame(r1, r2)
    colnames(observed) <- colnames(choices)
    assignments <- match(observed, choices)

    STATUS <- function(R1, R2) (R1!=0L) + 2*(R2!=0L)
    status <- STATUS(r1, r2)

    # Handling the inverted case.
    extra.qc <- list()
    if (randomized) {
        r1.alt <- .solo_identifier(fastq[[2]], choices=condensed[[1]], flank5=flank5[1], flank3=flank3[1],
            substitutions=substitutions[1], insertions=insertions[1], deletions=deletions[1], total.edits=total.edits[1],
            strand=strand[1])
        r2.alt <- .solo_identifier(fastq[[1]], choices=condensed[[2]], flank5=flank5[2], flank3=flank3[2],
            substitutions=substitutions[2], insertions=insertions[2], deletions=deletions[2], total.edits=total.edits[2],
            strand=strand[2])

        observed.alt <- DataFrame(r1.alt, r2.alt)
        colnames(observed.alt) <- colnames(choices)
        assignments.alt <- match(observed.alt, choices)

        old.status <- status
        new.status <- STATUS(r1.alt, r2.alt)
        status <- pmax(old.status, new.status)

        # Ignoring ambiguous orientations where both orientations yield a valid match 
        # (that is not the same match, e.g., in the case of symmetric constructs!)
        has.both <- !is.na(assignments) & !is.na(assignments.alt) & assignments!=assignments.alt
        replace <- is.na(assignments) & !is.na(assignments.alt)
        assignments[replace] <- assignments.alt[replace]
        assignments[has.both] <- NA

        extra.qc <- list(
             original.orientation=sum(!is.na(assignments) & is.na(assignments.alt)),
             reverse.orientation=sum(replace),
             ambiguous.orientation=sum(has.both)
        )

        if (include.invalid) {
            # More generous replacement than just 'replace', as we ensure 
            # that invalid pairs are swapped in.
            replace2 <- replace | old.status!=3L & new.status==3L
            observed[replace2,] <- observed.alt[replace2,]
        }
    }

    output <- original
    output$counts <- tabulate(assignments, nbins=nrow(original))

    invalid.pair <- is.na(assignments) & status==3L
    if (include.invalid) {
        # Adding all our invalid observed friends.
        invalids <- observed[invalid.pair,]
        m <- selfmatch(invalids)
        keep <- !duplicated(m)
        invalids <- invalids[keep,,drop=FALSE]
        invalids[,1] <- original[invalids[,1],1]
        invalids[,2] <- original[invalids[,2],2]
        invalids$counts <- countMatches(m[keep], m)
        invalids$valid <- rep(FALSE, nrow(invalids))
        output$valid <- rep(TRUE, nrow(output))
        output <- rbind(output, invalids)
    }

    metadata(output) <- c(
        list(
            none=sum(status==0L), 
            barcode1.only=sum(status==1L),
            barcode2.only=sum(status==2L), 
            invalid.pair=sum(invalid.pair)
        ),
        extra.qc
    )

    output
}

#' @importFrom ShortRead FastqStreamer yield sread
.solo_identifier <- function(fastq, choices, flank5, flank3, substitutions, insertions, deletions, total.edits, strand) {
    strand <- match.arg(strand, c("original", "both", "reverse"))
    use.forward <- strand %in% c("original", "both")
    use.reverse <- strand %in% c("reverse", "both")

    ptr <- setup_barcodes_single(c(flank5, flank3), choices, substitutions, insertions, deletions, total.edits)
    incoming <- FastqStreamer(fastq) 
    on.exit(close(incoming))

    collected <- list()
    counter <- 1L
    while (length(fq <- yield(incoming))) {
        collected[[counter]] <- identify_barcodes_single(sread(fq), ptr, use.forward, use.reverse)
        counter <- counter + 1L
    }

    unlist(collected)
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
            paths2=vapply(files, "[", i=2, ""), 
            output, nmapped=colSums(mat)))

    if (withDimnames) {
        colnames(se) <- basename(se$paths1)
    }
    se 
}
