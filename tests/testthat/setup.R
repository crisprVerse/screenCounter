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
    collected
}


