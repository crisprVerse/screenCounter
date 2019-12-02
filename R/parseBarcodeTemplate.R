#' Parse barcode template
#' 
#' Parse a barcode template to identify variable regions based on the run of N's.
#' 
#' @param template String containing template sequence of a barcode.
#' Variable regions should be marked with N's.
#' 
#' @return A list containing:
#' \itemize{
#' \item \code{variable}, a \linkS4class{DataFrame} containing the position and length of each run of N's.
#' \item \code{constant}, a character vector of constant regions flanking and separating the variable regions.
#' }
#'
#' @details
#' The barcode template should contain runs of N's to mark the variable regions.
#' The first run of N's is the first variable region, the second run of N's is the second variable region, and so on.
#' The template is \dQuote{realized} into a barcode when the N's are replaced with actual DNA sequence.
#' The use of a template provides a convenient format to express the general structure of the barcode while avoiding confusion about barcode-specific variable regions.
#' 
#' @author Aaron Lun
#' 
#' @examples
#' # Single spacer:
#' parseBarcodeTemplate("AAAANNNNNNNGGGG")
#'
#' # Double spacer:
#' parseBarcodeTemplate("AAAANNNNCCCCNNGGGG")
#' @export
parseBarcodeTemplate <- function(template) {
    positions <- .split_template(template)
    n.pos <- positions$pos
    n.len <- positions$len

    # Extracting constant regions.
    last.pos <- 1L
    constants <- character(length(n.pos)+1L)
    for (i in seq_along(n.pos)) {
        constants[i] <- substring(template, last.pos, n.pos[i]-1)
        last.pos <- n.pos[i] + n.len[i]
    }
    constants[length(n.pos)+1] <- substring(template, last.pos, nchar(template))

    list(variable=DataFrame(pos=n.pos, len=n.len), constant=constants)
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



