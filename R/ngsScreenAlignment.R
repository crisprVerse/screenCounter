#' @importFrom stringr str_extract
.getFastqFiles <- function(ngs, runs=NULL){
    fastq.df <- ngsprojFastq(ngs,collapse=FALSE)
    fastq.df$run <- str_extract(fastq.df$READ1_FILE,"R[0-9]+")
    if (!is.null(runs)){
        if (!all(runs %in% fastq.df$run)){
            availRuns <- paste0(unique(fastq.df$run), collapse=" ")
            stop(paste0("Only runs ", availRuns, " are available."))
        } else {
            fastq.df <- fastq.df[fastq.df$run %in% runs,,drop=FALSE]
        }
    }
    if (sum(is.na(fastq.df$READ2_FILE))>0){
        fastqs <- fastq.df$READ1_FILE
    } else {
        fastqs <- c(fastq.df$READ1_FILE, fastq.df$READ2_FILE)
    }
    fastqs <- gsub("gne/research", "gne", fastqs)
    fastqs 
}


#' Process reads from CRISPR and ORF screens
#' 
#' Process reads from CRISPR and ORF screens and return a SummarizedExperiment object. 
#'
#' @param ngs String specyfing the NGS number (eg. "ngs3096").
#' @param libname String specifying the name of the CRISPR/ORF library. Must be available in crisprAnnotation.
#' @param flank5 String containing the constant sequence on the 5' flank of the variable region.
#' @param flank3 String containing the constant sequence on the 3' flank of the variable region.
#' @param runs Character string indicating which sequencing runs should be considered.
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
#'     flank5="TACCG", flank3="GTTT"
#' )
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
ngsScreenAlignment <- function(ngs, libname, flank5, flank3, runs=NULL,...){

    files <- .getFastqFiles(ngs, runs=runs)
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




#' Process reads for cell barcoding experiment
#' 
#' Process reads for cell barcoding experiment using lentiviral delivery.
#'
#' @param ngs String specyfing the NGS number (eg. "ngs3096").
#' @param choices A \linkS4class{List} of character vectors, one per variable region in \code{template}.
#' The first vector should contain the potential sequences for the first variable region, 
#' the second vector for the second variable region and so on.
#' @param template A template for the barcode structure, see \code{?\link{parseBarcodeTemplate}} for details.
#' @param runs Character string indicating which sequencing runs should be considered.
#' @param ... Further arguments to pass to \code{matrixOfSingleBarcodes}.
#'
#' @return 
#  Will return a \linkS4class{SummarizedExperiment} object containing:
#' \itemize{
#' \item An integer matrix named \code{"counts"}, where each column corresponds to a SAMID.
#' \item Row metadata containing the barcodes information
#' \item Column metadata containing sample information as annotated in HiTS-LIMS, as well as;
#' an integer vector \code{nreads}, containing the total number of reads in each file;
#' and \code{nmapped}, containing the number of reads assigned to a combo barcode in; 
#' the output count matrix.
#' }
#'
#' @author Jean-Philippe Fortin
#' @examples
#' # Processing a tumor barcoding experiment using Cellecta50M lib:
#' \dontrun{
#' library(crisprAnnotation)
#' barcodes <- getBarcodeLibrary("barcodes.cellecta.celltracker.50M")
#' varRegion <- paste0(rep("N",18), collapse="")
#' template <- paste0("GTTCG", 
#'    varRegion,"TTCG",varRegion, 
#'    "TTCGG"
#' )
#' se <- getBarcodeCounts(ngs, list(barcodes, barcodes), template)
#' }
#' 
#' @export
#' @importFrom gneDB ngsprojFastq
#' @importFrom gneDB annotateSAMIDs
#' @importFrom gp.sa.utils RosalindParam
#' @importFrom stringr str_extract
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment assays
ngsComboBarcodingAlignment <- function(ngs, choices, template, runs=runs,...){
    if (length(choices)!=2){
        stop("Only dual barcodes supported at the moment. ")
    }
    files <- .getFastqFiles(ngs,runs=runs)
    se <- matrixOfComboBarcodes(files,
        choices=choices, 
        template=template,
        strand="both",
        BPPARAM=RosalindParam(length(files)),
        withDimnames=TRUE,
        ... 
    )      
    counts <- assays(se)[[1]]

    # Sample key:
    sample <- colnames(counts)
    sample <- str_extract(sample, "SAM[0-9]+_R[1-2]")
    sample <- gsub("_R1|_R2","",sample)
    samples <- unique(sample)
    Y <- lapply(samples, function(x){
        wh <- which(sample==x)
        rowSums(counts[,wh, drop=FALSE], na.rm=TRUE)
    }) 
    Y <- do.call(cbind, Y)
    colnames(Y) <- samples
    
     # Creating summarized colData:
    analyzed <- as.matrix(colData(se)[, c("nreads", "nmapped")])
    Z <- lapply(samples, function(x){
        wh <- which(sample==x)
        colSums(analyzed[wh,, drop=FALSE], na.rm=TRUE)
    }) 
    Z <- do.call(rbind, Z)
    rownames(Z) <- samples
   
    # Let's get the phenotype information:
    sinfo <- annotateSAMIDs(rownames(Z))[colnames(Y),]
    sinfo[is.na(sinfo)] <- ""
    sinfo <- cbind(sinfo, Z)

    # Let's create the rowData:
    ann <- rowData(se)
    colnames(ann) <- c("barcode1", "barcode2")
    rownames(ann) <- paste0("Barcode_", 1:nrow(ann))
    rownames(Y) <- rownames(ann)
    # Creating new SE:
    se <- SummarizedExperiment(Y, 
        colData=sinfo, 
        rowData=ann, 
        metadata=list(ngs=ngs)
    )
    se <- se[,order(colnames(se))]
    se
}



