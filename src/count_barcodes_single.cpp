#include "Rcpp.h"

#include "hash_sequence.h"
#include "build_hash.h"

/* Single-end guide parser. */

// [[Rcpp::export(rng=false)]]
SEXP count_barcodes_single(SEXP seqs, SEXP guides) {
    auto hash=build_hash(guides);
    const auto& reference=hash.first;
    const auto& all_lens=hash.second;

    // Running through the sequences and matching it to the guides.
    Rcpp::StringVector Seqs(seqs);
    Rcpp::IntegerVector output(Seqs.size(), -1);

    for (size_t i=0; i<Seqs.size(); ++i) {
        Rcpp::String s=Seqs[i];
        const char* sptr=s.get_cstring();
        const size_t len=Rf_length(s.get_sexp());

        for (auto curlen : all_lens) {
            if (curlen > len) { break; }

            auto curstr=hash_sequence(sptr, curlen);
            auto it=curhash.find(curstr);
            if (it!=curhash.end()) {
                output[i]=it->second;
                break;
            }

            size_t j=0;
            while (j + curlen < len) {
                shift_sequence(curstr, curlen, sptr[j+curlen]);
                auto it=curhash.find(curstr);
                if (it!=curhash.end()) {
                    output[i]=it->second;
                    break;
                }
                ++j;
            }
        }
    }

    return output;
}
