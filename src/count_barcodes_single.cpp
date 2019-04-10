#include "Rcpp.h"

#include "hash_sequence.h"
#include "build_hash.h"

/* Single-end guide parser. */

typedef std::pair<seqhash, std::vector<int> > hashinfo;

// [[Rcpp::export(rng=false)]]
SEXP setup_barcodes_single(SEXP guides) {
    hashinfo* ptr = new hashinfo(build_hash(guides));
    return Rcpp::XPtr<hashinfo>(ptr, true);
}

// [[Rcpp::export(rng=false)]]
SEXP count_barcodes_single(SEXP seqs, SEXP xptr) {
    Rcpp::XPtr<hashinfo> ptr(xptr);
    const auto& reference=ptr->first;
    const auto& all_lens=ptr->second;

    // Running through the sequences and matching it to the guides.
    Rcpp::StringVector Seqs(seqs);
    Rcpp::IntegerVector output(Seqs.size());

    for (size_t i=0; i<Seqs.size(); ++i) {
        Rcpp::String s=Seqs[i];
        const char* sptr=s.get_cstring();
        const size_t len=Rf_length(s.get_sexp());

        for (auto curlen : all_lens) {
            if (curlen > len) { break; }

            auto curstr=hash_sequence(sptr, curlen);
            int nvalid=is_valid(sptr, curlen);
            if (nvalid==curlen) {
                auto it=reference.find(curstr);
                if (it!=reference.end()) {
                    output[i]=it->second;
                    break;
                }
            }

            size_t start=0, end=curlen;
            while (end < len) {
                nvalid-=is_valid(sptr[start]);
                shift_sequence(curstr, curlen, sptr[end]);
                nvalid+=is_valid(sptr[end]);

                if (nvalid==curlen) {
                    auto it=reference.find(curstr);
                    if (it!=reference.end()) {
                        output[i]=it->second+1; // Get back to 1-indexed.
                        break;
                    }
                }

                ++end;
                ++start;
            }
        }
    }

    return output;
}
