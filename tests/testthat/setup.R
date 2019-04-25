# Basic hashing utilities to be used in all other tests.

BASES <- c("A", "C", "G", "T")

GENERATE_RANDOM_SEQ <- function(n) {
    paste(sample(BASES, n, replace=TRUE), collapse="")
}

MANUAL_HASH <- function(barcode, split=FALSE) {
    if (split) {
        B <- barcode
    } else { 
        B <- strsplit(barcode, "")[[1]]
    }

    collected <- numeric(0)
    while (length(B)) {
        word <- head(B, 16)
        as_num <- c(A=0, C=1, G=2, T=3)[word] * 4^(seq_along(word)-1)
        collected <- c(collected, sum(as_num))
        B <- tail(B, -16)
    }

    # Adding on the missing half of the 64-bit integer.
    if (length(collected)%%2==1) {
        collected <- c(collected, 0)
    }
    collected
}

ADD_FLANKS <- function(barcodes, fname, nleft=50, nright=50) {
    N <- length(barcodes)
    left <- vapply(sample(nleft, N, replace=TRUE), GENERATE_RANDOM_SEQ, FUN.VALUE="")
    right <- vapply(sample(nright, N, replace=TRUE), GENERATE_RANDOM_SEQ, FUN.VALUE="")
    barcodes2 <- paste0(left, barcodes, right)
    names(barcodes2) <- seq_len(N)
    writeXStringSet(DNAStringSet(barcodes2), filepath=fname, format="fastq")
}
