#include "Rcpp.h"
#include "utils.h"
#include <stdexcept>

kaori::BarcodePool format_pointers(const Rcpp::CharacterVector& options) {
    size_t size = 0, num_opts = options.size();

    std::vector<const char*> ptrs(num_opts);
    for (size_t o = 0; o < num_opts; ++o) {
        auto current = Rcpp::String(options[o]);

        size_t curlen = LENGTH(current.get_sexp());
        if (!o) {
            size = curlen;
        } else if (size != curlen) {
            throw std::runtime_error("variable regions should all have the same length (" + std::to_string(size) + ")");
        }

        ptrs[o] = current.get_cstring();
    }

    return kaori::BarcodePool(std::move(ptrs), size);
}

kaori::SearchStrand to_strand(bool reverse) {
    if (reverse) {
        return kaori::SearchStrand::REVERSE;
    } else {
        return kaori::SearchStrand::FORWARD;
    }
}

kaori::SearchStrand to_strand(int strand) {
    if (strand == 0) {
        return kaori::SearchStrand::FORWARD;
    } else if (strand == 1) {
        return kaori::SearchStrand::REVERSE;
    } else {
        return kaori::SearchStrand::BOTH;
    }
}
