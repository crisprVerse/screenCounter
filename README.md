# Counting barcodes in functional screens

|Environment|Status|
|---|---|
|[BioC-release](https://bioconductor.org/packages/release/bioc/html/screenCounter.html)|[![Release OK](https://bioconductor.org/shields/build/release/bioc/screenCounter.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/screenCounter/)|
|[BioC-devel](https://bioconductor.org/packages/devel/bioc/html/screenCounter.html)|[![Devel OK](https://bioconductor.org/shields/build/devel/bioc/screenCounter.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/screenCounter/)|

**screenCounter** quantifies the frequency of barcodes in sequencing data from functional screens,
yielding a count matrix that can be subjected to further statistical analyses, e.g., with [**edgeR**](https://bioconductor.org/packages/edgeR) and contemporaries.
We support single barcode, combinatorial barcode and dual barcode designs, where each construct is assumed to have constant sequence(s) flanking one or more variable regions of known length.
The package also provides options for mismatches, strand-specific searching and parallelization via [**BiocParallel**](https://bioconductor.org/packages/BiocParallel).

_For users:_ check out the documentation at the [Bioconductor landing page](https://bioconductor.org/packages/screenCounter) for installation and usage instructions.

_For developers:_ this R package is just a wrapper around the [**kaori**](https://github.com/crisprVerse/kaori) C++ library that does the heavy lifting.
Pull requests for new barcode designs should be made in **kaori** first, and then propagated back down to **screenCounter**.
