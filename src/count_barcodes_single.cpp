#include "Rcpp.h"

#include "hash_string.h"
#include <unordered_map>
#include <deque>

/* Single-end guide parser. */

// [[Rcpp::export(rng=false)]]
SEXP count_barcodes_single(SEXP seqs, SEXP guides) {
    std::vector<int> len2hash;
    std::deque<std::unordered_map<std::u32string, int> > all_hashes;

    // Constructing the guide hash map.
    Rcpp::StringVector Guides(guides);
    for (size_t i=0; i<Guides.size(); ++i) {
        Rcpp::String g=Guides[i];

        const size_t len=Rf_length(g.get_sexp());
        if (len >= len2hash.size()) {
            len2hash.resize(len+1, -1);
        }

        if (len2hash[len]==-1) {
            len2hash[len]=all_hashes.size();
            all_hashes.resize(all_hashes.size()+1);
        }

        all_hashes[len2hash[len]][hash_sequence(g.get_cstring(), len)] = i;
    }

    std::deque<std::pair<int, int> > all_lens;
    for (size_t i=0; i<len2hash.size(); ++i) {
        if (len2hash[i]!=-1) {
            all_lens.push_back(std::make_pair(i, len2hash[i]));
        }
    }

    // Running through the sequences and matching it to the guides.
    Rcpp::StringVector Seqs(seqs);
    Rcpp::IntegerVector output(Seqs.size(), -1);

    for (size_t i=0; i<Seqs.size(); ++i) {
        Rcpp::String s=Seqs[i];
        const char* sptr=s.get_cstring();
        const size_t len=Rf_length(s.get_sexp());

        for (auto current : all_lens) {
            if (current.first > len) { break; }
            size_t curlen=current.first;
            const auto& curhash=all_hashes[current.second];

            auto curstr=hash_sequence(sptr, curlen);
            auto it=curhash.find(curstr);
            if (it!=curhash.end()) {
                output[i]=it->second;
                break; // TODO: keep looking if not a perfect match.
            }

            size_t j=0;
            while (j + curlen < len) {
                shift_string(curstr, curlen, sptr[j+curlen]);
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
