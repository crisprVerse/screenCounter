# Counting barcodes in functional screens

Pretty much what it says on the tin.
The **screenCounter** package provides functions to count barcodes in sequencing data from functional screens,
yielding a count matrix that can be subjected to further statistical analyses, e.g., with [**edgeR**](https://bioconductor.org/packages/edgeR) and contemporaries.
We support single barcode, combinatorial barcode and dual barcode designs, where each construct is assumed to have constant sequence(s) flanking one or more variable regions.
The package also provides options for mismatches, strand-specific searching and parallelization via [**BiocParallel**](https://bioconductor.org/packages/screenCounter).
