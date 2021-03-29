#include "utils.h"
#include <algorithm>

void reverse_complement (char* seq, size_t len) {
    for (size_t i=0; i<len; ++i) {
        char& current=seq[i];
        switch(current) {
            case 'A': case 'a':
                current='T'; break;
            case 'C': case 'c':
                current='G'; break;
            case 'G': case 'g':
                current='C'; break;
            case 'T': case 't':
                current='A'; break;
        }
    }
    std::reverse(seq, seq+len);
    return;
}

/* Utility functions to check that a string or character is valid,
 * i.e., contains only ACTG (or lower case) bases.
 */

bool is_valid(char base) {
    switch (base) {
        case 'A': case 'a': case 'C': case 'c': case 'G': case 'g': case 'T': case 't':
            return true;
    };
    return false;
}

int is_valid(const char* ptr, size_t n) {
    int valid=0;
    for (size_t i=0; i<n; ++i) {
        valid+=is_valid(ptr[i]);
    }
    return valid;
}

size_t check_length (const Rcpp::StringVector& incoming) {
    size_t reflen = -1;
    for (size_t v = 0; v < incoming.size(); ++v) {
        Rcpp::String vc = incoming[v];
        const size_t vlen = Rf_length(vc.get_sexp());

        if (reflen == -1) {
            reflen = vlen;
        } else if (reflen != vlen) {
            throw std::runtime_error("barcodes are of variable length");
        }
    }
    if (reflen < 0) {
        throw std::runtime_error("no barcodes supplied");
    }
    return reflen;
}
