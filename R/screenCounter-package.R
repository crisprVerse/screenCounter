#' The screenCounter package
#'
#' This package provides methods to counting barcodes from read sequences in high-throughput sequencing screen data sets.
#' It does so by loading in chunks of sequences from FASTQ files using \pkg{ShortRead}
#' and then matching the barcode template to each sequence using a rolling hash (implemented in C++, inspired by Colin Watanabe's code).
#' This process is performed across several files using a range of parallelization schemes available in \pkg{BiocParallel}.
#' We return the resulting count matrix and any feature annotations in a \linkS4class{SummarizedExperiment} object.
#' Currently, single barcodes (in single- or paired-end data) and combinatorial barcodes are supported.
#'
#' @author Aaron Lun
#'
#' @name screenCounter-package
#' @importFrom Rcpp sourceCpp
#' @importFrom IRanges IRanges
#' @importFrom Biostrings DNAStringSet
#' @importFrom XVector XRaw
#' @useDynLib screenCounter
NULL
