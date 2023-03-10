#include "Rcpp.h"
#include "utils.h"
#include <stdexcept>

kaori::BarcodePool format_pointers(const Rcpp::CharacterVector& options) {
    size_t size = 0;

    std::vector<const char*> ptrs(options.size());
    for (size_t o = 0; o < options.size(); ++o) {
        auto current = Rcpp::String(options[o]);

        if (o) {
            if (size != LENGTH(current.get_sexp())) {
                throw std::runtime_error("variable regions should all have the same length (" + std::to_string(size) + ")");
            }
        } else {
            size = LENGTH(current.get_sexp());
        }

        ptrs[o] = current.get_cstring();
    }

    return kaori::BarcodePool(std::move(ptrs), size);
}
