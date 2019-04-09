#' Create barcodes from a template
#' 
#' Create barcodes from a template template by filling in N's with possible sequences.
#' 
#' @param template String containing template sequence of a barcode.
#' Variable regions should be marked with N's.
#' @param choices A \linkS4class{DataFrame} of potential choices for the variable regions.
#' Each row represents a distinct barcode, and each column represents one of the variable regions
#' (in their order on \code{template}, from 5' to 3').
#' 
#' @return Character vector containing realized barcodes.
#' This is named with the row names of \code{choices}.
#'
#' @details 
#' The choices for the variable regions do not have to have the same length as in \code{template}.
#' This potentially allows for variable-length barcodes, if such a thing is desirable.
#'
#' @author Aaron Lun
#' 
#' @examples
#' # Single spacer:
#' createBarcodes("AAAANNNNCCCCNNGGGG", 
#'    DataFrame(c("T", "TT", "TTT")))
#'
#' # Double spacer:
#' createBarcodes("AAAANNNNCCCCNNGGGG", 
#'    DataFrame(c("T", "TT", "TTT"), c("TTTT","TT", "T")))
#' @export
createBarcodes <- function(template, choices) {
    if (nrow(choices)==0L) {
        return(character())
    }

    positions <- .split_template(template)
    n.pos <- positions$pos
    n.len <- positions$len

    # Validating columns.
    nvariables <- length(n.pos)
    if (nvariables!=ncol(choices)) {
        stop("'ncol(choices)' is not equal to the number of stretches of N's")
    }

    # Interleaving constant and variable regions.
    collected <- vector("list", 1L+2L*nvariables)
    last.pos <- 1L
    for (i in seq_along(n.pos)) {
        collected[[i*2L-1L]] <- substring(template, last.pos, n.pos[i]-1)
        collected[[i*2L]] <- choices[,i]
        last.pos <- n.pos[i] + n.len[i]
    }

    collected[[length(n.pos)*2L + 1L]] <- substring(template, last.pos, nchar(template))

    output <- do.call(paste0, collected)
    names(output) <- rownames(choices)
    output
}

.split_template <- function(template) {
    has.n <- gregexpr("N+", template)[[1]]
    n.pos <- as.integer(has.n)
    if (identical(n.pos, -1L)) {
        stop("no stretches of Ns detected in 'template'")
    } 

    n.len <- attr(has.n, "match.length")
    list(pos=n.pos, len=n.len)
}



