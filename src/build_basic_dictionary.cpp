#include "build_basic_dictionary.h"

basic_dictionary build_basic_dictionary (Rcpp::StringVector guides) {
    basic_dictionary mapping;
    for (size_t i=0; i<guides.size(); ++i) {
        Rcpp::String g=guides[i];
        const size_t len=Rf_length(g.get_sexp());
        const char* ptr=g.get_cstring();
        auto curhash=hash_sequence(ptr, len);

        auto it=mapping.find(curhash);
        if (it!=mapping.end()) {
            throw std::runtime_error("duplicated barcode sequences");
        } else {
            mapping[curhash]=i;
        }
    }
    return mapping;
}
