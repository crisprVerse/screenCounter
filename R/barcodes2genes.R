#' Barcode to gene statistics
#'
#' Combine barcode-level differential abundance statistics into gene-level statistics.
#' 
#' @param x A \linkS4class{DataFrame} of per-barcode results. 
#' @param genes A character vector or factor containing the gene ID for each barcode.
#' @param pval.col String specifying the column of \code{x} containing the p-values.
#' @param lcpm.col String specifying the column of \code{x} containing the log-CPMs.
#' @param lfc.regex String containing a regular expression to match the column names of \code{x} containing the log-fold changes.
#'
#' @return
#' A \linkS4class{DataFrame} containing per-gene statistics, with one row per gene.
#' 
#' @details
#' Statistics from all barcodes of a particular gene are defined as follows:
#' \itemize{
#' \item \code{NBarcodes}: the number of barcodes for a given gene.
#' \item \code{NLogFC.up}, \code{NLogFC.down}: the number of barcodes with log2-fold changes greater than 0.5 or less than -0.5, respectively.
#' \item \code{PValue}: the p-value for each gene.
#' This is computed by combining all p-values from the relevant barcodes with Simes' method.
#' Here, the joint null hypothesis is that no barcodes are differentially abundant for a gene.
#' \item \code{FDR}: the Benjamini-Hochberg adjusted p-value for each gene.
#' \item \code{Direction}: the overall direction of the change in abundance across barcodes for each gene.
#' This is determined from the direction of the tests with p-values small enough to contribute to the \code{PValue}.
#' See \code{\link{combineTests}} for more details.
#' \item \code{LogFC}: the log-fold change of the most significant barcode (i.e., lowest p-value) for each gene.
#' \item \code{LogCPM}: the log-CPM of the most significant barcode (i.e., lowest p-value) for each gene.
#' }
#' 
#' If ANOVA-like contrasts are supplied in \code{x}, \code{Direction} is not reported, 
#' and each \code{LogFC} column name will be suffixed by the name of the relevant column of the contrast.
#' 
#' @author Aaron Lun
#' @seealso
#' \code{\link{combineTests}} and \code{\link{getBestTest}}, for the actual statistical calculations.
#'
#' @examples
#' example(DAScreenStatFrame, echo=FALSE)
#'
#' genes <- sample(LETTERS[1:3], nrow(Y), replace=TRUE)
#' barcodes2genes(Y, genes=genes)
#' 
#' @export
#' @importFrom csaw combineTests getBestTest
#' @importFrom S4Vectors DataFrame
barcodes2genes <- function(x, genes, pval.col="PValue", lcpm.col="LogCPM", lfc.regex="^LogFC") {
    lost <- is.na(genes)
    all.genes <- sort(unique(genes[!lost]))
    keep <- !is.na(x[,pval.col]) & !lost
    genes <- genes[keep]
    x <- x[keep,,drop=FALSE]

    per.gene <- combineTests(genes, x, pval.col=pval.col, fc.col=grep(lfc.regex, colnames(x)))
    colnames(per.gene)[1] <- 'NBarcodes'
    colnames(per.gene)[colnames(per.gene)=="direction"] <- "Direction"
    is.lfc <- grep(lfc.regex, colnames(per.gene))
    colnames(per.gene)[is.lfc] <- paste0("N", colnames(per.gene)[is.lfc])

    best <- getBestTest(genes, x, pval.col=pval.col)
    best <- best[,c(which(colnames(best)==lcpm.col), grep(lfc.regex, colnames(best))),drop=FALSE]

    output <- cbind(per.gene, best)
    output <- output[match(all.genes, rownames(output)),]
    rownames(output) <- all.genes
    
    DataFrame(output)
}
