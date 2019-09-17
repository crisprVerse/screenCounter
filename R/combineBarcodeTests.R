#' Combine barcode statistics
#'
#' Combine barcode-level differential abundance statistics into gene-level statistics.
#' 
#' @param x A \linkS4class{DataFrame} of per-barcode results. 
#' @param genes A character vector or factor containing the gene ID for each barcode.
#' @param pval.col String specifying the column of \code{x} containing the p-values.
#' @param lcpm.col String specifying the column of \code{x} containing the log-CPMs.
#' @param lfc.regex String containing a regular expression to match the column names of \code{x} containing the log-fold changes.
#' @param method String specifying the method to use to combine barcode-level p-values into a per-gene p-value.
#' @param min.sig.n Integer scalar containing the minimum number of significant barcodes when \code{method="holm-min"}.
#' @param min.sig.prop Numeric scalar containing the minimum proportion of significant barcodes when \code{method="holm-min"}.
#'
#' @return
#' A \linkS4class{DataFrame} containing per-gene statistics, with one row per gene and the following fields:
#' \itemize{
#' \item \code{NBarcodes}: the number of barcodes for a given gene.
#' \item \code{NLogFC.up}, \code{NLogFC.down}: the number of barcodes with log2-fold changes greater than 0.5 or less than -0.5, respectively.
#' \item \code{PValue}: the p-value for each gene, combined from the per-barcode p-values.
#' \item \code{FDR}: the Benjamini-Hochberg adjusted p-value for each gene.
#' \item \code{Direction}: the overall direction of the change in abundance across barcodes for each gene.
#' This is determined from the direction of the tests with p-values small enough to contribute to the \code{PValue}.
#' See \code{\link{combineTests}} for more details.
#' \item \code{LogFC}: the log-fold change of the most significant barcode (i.e., lowest p-value) for each gene.
#' \item \code{LogCPM}: the log-CPM of the most significant barcode (i.e., lowest p-value) for each gene.
#' }
#' 
#' @details
#' When \code{method="simes"}, the per-gene p-value is computed by combining all p-values from the relevant barcodes with Simes' method.
#' Here, the global null hypothesis is that no barcodes are differentially abundant for a gene.
#' With this approach, a gene can obtain a low p-value even if only a few guides (or just one guide) are strongly differentially abundant.
#' This is the most sensitive approach but is also more susceptible to off-target effects.
#' It is implemented with the \code{\link{combineTests}} function from the \pkg{csaw} package.
#'
#' When \code{method="holm"}, we apply a Holm correction across all barcodes associated with a given gene,
#' and then set the combined p-value to the \eqn{k}-th largest Holm-corrected barcode-level p-value.
#' Here, the global null hypothesis is that fewer than \eqn{k} barcodes are differentially abundant for a gene.
#' Thus, a gene it can only obtain a low p-value if it has at least \eqn{k} differentially abundant barcodes.
#' 
#' To define \eqn{k} for any gene, we take the larger of \code{min.sig.n} and the product of \code{min.sig.prop} with the number of barcodes for that gene.
#' This is then capped at the number of barcodes for that gene.
#' By setting \eqn{k > 1}, we require some level of agreement between barcodes and reduce the influence of one off-target guide on the gene-level inference.
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
#' combineBarcodeTests(Y, genes=genes)
#' 
#' @seealso
#' \code{\link{barcodeSetTest}}, for another way of consolidating barcode statistics to gene-level results.
#'
#' @export
#' @importFrom csaw combineTests getBestTest
#' @importFrom S4Vectors DataFrame
#' @importFrom stats p.adjust
combineBarcodeTests <- function(x, genes, pval.col="PValue", lcpm.col="LogCPM", lfc.regex="^LogFC", 
    method=c("simes", "holm-min"), min.sig.n=3, min.sig.prop=0.4)
{
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

    if (match.arg(method)!="simes") {
        # Computing the Holm minimal p-value.
        barcode.p <- split(x[,pval.col], genes)
        for (i in names(barcode.p)) {
            p <- p.adjust(barcode.p[[i]], method="holm")    
            barcode.p[[i]] <- sort(p)[min(length(p), max(min.sig.n, ceiling(min.sig.prop * length(p))))]
        }
        per.gene[[pval.col]] <- unlist(barcode.p)[rownames(per.gene)]
        per.gene$FDR <- p.adjust(per.gene[[pval.col]], method="BH")
    }

    best <- getBestTest(genes, x, pval.col=pval.col)
    best <- best[,c(which(colnames(best)==lcpm.col), grep(lfc.regex, colnames(best))),drop=FALSE]
    output <- cbind(per.gene, best)

    output <- output[match(all.genes, rownames(output)),]
    rownames(output) <- all.genes
    DataFrame(output)
}
