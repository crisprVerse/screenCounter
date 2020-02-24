#' Process reads from CRISPR and ORF screens
#' 
#' Process reads from CRISPR and ORF screens and return a SummarizedExperiment object. 
#'
#' @param ngs String specyfing the NGS number (eg. "ngs3096").
#' @param libname String specifying the name of the CRISPR/ORF library. Must be available in crisprAnnotation.
#' @param flank5 String containing the constant sequence on the 5' flank of the variable region.
#' @param flank3 String containing the constant sequence on the 3' flank of the variable region.
#' @param ... Further arguments to pass to \code{matrixOfSingleBarcodes}.
#'
#' @return 
#  Will return a \linkS4class{SummarizedExperiment} object containing:
#' \itemize{
#' \item An integer matrix named \code{"counts"}, where each column corresponds to a SAMID.
#' \item Row metadata containing the CRISPR library annotation.
#' \item Column metadata containing sample information as annotated in HiTS-LIMS, as well as;
#' an integer vector \code{nreads}, containing the total number of reads in each file;
#' and \code{nmapped}, containing the number of reads assigned to a barcode in the output count matrix.
#' }
#'
#' @author Jean-Philippe Fortin
#' @examples
#' # Processing CRISPRa screen NGS3096 (Pilot screen with Transylvania library)
#' \dontrun{
#' ngsScreenAlignment(ngs="ngs3096",
#'     libname="crispra.cas9.human.transylvania1",
#'     flank5="TACCG", flank3="GTTT")
#' }
#' @export
#' @importFrom gneDB ngsprojFastq
#' @importFrom gneDB annotateSAMIDs
#' @importFrom crisprAnnotation listCrisprLibraries
#' @importFrom crisprAnnotation getCrisprLibraryAnnotation
#' @importFrom gp.sa.utils RosalindParam
#' @importFrom stringr str_extract
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment assays
ngsScreenAlignment <- function(ngs, libname, flank5, flank3, ...){

	.getFastqFiles <- function(ngs){
	    fastq.df <- ngsprojFastq(ngs,collapse=FALSE)
	    if (sum(is.na(fastq.df$READ2_FILE))>0){
	        fastqs <- fastq.df$READ1_FILE
	    } else {
	        fastqs <- c(fastq.df$READ1_FILE, fastq.df$READ2_FILE)
	    }
	    fastqs <- gsub("gne/research", "gne", fastqs)
	    fastqs 

	}
	files <- .getFastqFiles(ngs)
	if (!libname %in% listCrisprLibraries()){
		stop("libname not found in crisprAnnotation")
	} 
	ann <- getCrisprLibraryAnnotation(libname)
	choices <- unique(ann$barcode)

  	se <- matrixOfSingleBarcodes(files,
        choices=choices, 
        flank5=flank5, 
        flank3=flank3,
        withDimnames=TRUE,
        BPPARAM=RosalindParam(length(files)),
        ...
    )
    counts <- assays(se)[[1]]

    #Merging reads1 with reads2:
    sample <- str_extract(colnames(counts), "SAM[0-9]+_R[1-2]")
    sample <- gsub("_R1|_R2","",sample)
    samples <- unique(sample)
    Y <- lapply(samples, function(x){
        wh <- which(sample==x)
        rowSums(counts[,wh, drop=FALSE], na.rm=TRUE)
    }) 
    Y <- do.call(cbind, Y)
    colnames(Y) <- samples

    #Subsetting annotation:
    ann <- ann[match(rownames(Y), ann$barcode),]
    rownames(Y) <- ann$id

    # Creating summarized colData:
    analyzed <- as.matrix(colData(se)[, c("nreads", "nmapped")])
    Z <- lapply(samples, function(x){
        wh <- which(sample==x)
        colSums(analyzed[wh,, drop=FALSE], na.rm=TRUE)
    }) 
    Z <- do.call(rbind, Z)
    rownames(Z) <- samples
   
    # Let's get the phenotype information:
    sinfo <- annotateSAMIDs(colnames(Y))[colnames(Y),]
    sinfo[is.na(sinfo)] <- ""
    sinfo <- cbind(sinfo, Z)
    
    # Creating new SE:
    se <- SummarizedExperiment(Y, 
        colData=sinfo, 
        rowData=ann, 
        metadata=list(ngs=ngs, lib=libname)
    )
    se <- se[,order(colnames(se))]
    se
}



