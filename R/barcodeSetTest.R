#' Perform a barcode set test
#'
#' Uses the \code{\link{fry}} function to perform a self-contained barcode set test for each gene.
#'
#' @param v Any object that can be passed to \code{\link{fry}}, e.g., a matrix of log-expression values, or an output of \code{\link{voom}}.
#' @param gene.ids A character vector specifying the gene targeted by each barcode.
#' @param all.genes A character vector specifying the universe of all genes.
#' @param ... Further arguments to pass to \code{\link{fry}}.
#' @param stats A data frame or \linkS4class{DataFrame} of per-barcode differential testing statistics, similar to that from \code{\link{topTable}}.
#' If provided, this should correspond to the entries in \code{gene.ids}.
#' @param lcpm.col String specifying the column of \code{stats} containing log-CPMs.
#' @param lfc.col String specifying the column of \code{stats} containing log-FCs.
#' @param subset Integer or logical vector specifying the subset of barcodes used to create \code{v}.
#'
#' @details
#' This function uses the \code{\link{fry}} machinery developed for self-contained gene set tests, and applies it to the set of barcodes for each gene.
#' The aim is to determine whether there is a consistent change in abundance across conditions for all barcodes associated with a given gene.
#' This contrasts with \code{\link{combineBarcodeTests}}, which is more willing to consider genes that only exhibit a change in abundance for barcode.
#' 
#' If \code{stats} is specified, some barcode-level statistics are added to the output DataFrame.
#' We report the median log-fold change and median log-CPM across all barcodes for each gene, using \code{lcpm.col} and \code{lfc.col}.
#' This aims to provide some additional information about the magnitude and abundance of the barcodes associated with each gene.
#' 
#' If \code{subset} is specified, rows of \code{v} are assumed to map to \code{gene.ids[subset]} and \code{res[subset,]}.
#' When combined with \code{all.genes}, this is useful for ensuring that the output DataFrame has one row per gene, for consistency in downstream reporting.
#' If a genes has no barcodes in the \code{subset}, its statistics are set to \code{NA}.
#' 
#' @return A \linkS4class{DataFrame} containing barcode set test statistics for each gene (row).
#' Columns include \code{"PValue"}, the p-value from the barcode set test; \code{"FDR"}, the adjusted p-value;
#' \code{"NBarcodes"}, the number of barcodes in \code{v} associated with a given gene; and
#' \code{"Direction"}, an overall direction of the change in abundance for the barcodes associated with a gene.
#' 
#' @author Aaron Lun
#' @examples
#' N <- 1000
#' mu <- 2^runif(N, 2, 10)
#' y <- matrix(rnbinom(N * 6, mu=mu, size=10), ncol=6)
#' g <- gl(2, 3)
#' ids <- sample(LETTERS, N, replace=TRUE)
#' 
#' library(limma)
#' design <- model.matrix(~g)
#' v <- voom(y, design=design)
#' 
#' output <- barcodeSetTest(v, design=design, contrast=2,
#'     gene.ids=ids)
#' output
#'
#' @seealso
#' \code{\link{combineBarcodeTests}}, another method for consolidating per-barcode statistics into per-gene results.
#'
#' @export
#' @importFrom stats median
#' @importFrom limma fry
#' @importFrom S4Vectors DataFrame
barcodeSetTest <- function(v, gene.ids, ..., 
    stats=NULL, lcpm.col="LogCPM", lfc.col="LogFC", subset=NULL,
    all.genes=sort(unique(gene.ids)))
{
    if (!is.null(subset)) {
        g <- gene.ids[subset]
    } else {
        g <- gene.ids
    }
    barcode.sets <- split(seq_along(g), g)
    gres <- fry(v, barcode.sets, ...) 

    if (!is.null(stats)) {
        if (!is.null(subset)) {
            stats <- stats[subset,]
        }
        gres[[lcpm.col]] <- vapply(barcode.sets, FUN=function(i) median(stats[i,lcpm.col], na.rm=TRUE), 0)
        gres[[lfc.col]] <- vapply(barcode.sets, FUN=function(i) median(stats[i,lfc.col], na.rm=TRUE), 0)
    }

    gres <- gres[match(all.genes, rownames(gres)),]
    rownames(gres) <- all.genes
    DataFrame(gres)
}
