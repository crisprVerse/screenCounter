#include "build_hash.h"
#include "hash_sequence.h"

/* Builds a hash index, returning it alongside a vector of lengths.
 * The latter is necessary to accommodate variable length barcodes.
 * It is okay to assume that there are no N's here, as any N's 
 * should have been eliminated by barcode construction at the R level.
 */

std::pair<seqhash, std::vector<int> > build_hash (Rcpp::StringVector guides) {
    seqhash reference;
    std::deque<bool> valid_lengths;

    for (size_t i=0; i<guides.size(); ++i) {
        Rcpp::String g=guides[i];
        const size_t len=Rf_length(g.get_sexp());

        // Checking that the length exists.
        if (len >= valid_lengths.size()) {
            valid_lengths.resize(len+1, false);
        }
        valid_lengths[len]=true;

        // Adding the true sequence, which is a high priority.
        const char* ptr=g.get_cstring();
        auto curhash=hash_sequence(ptr, len);

        auto it=reference.find(curhash);
        if (it!=reference.end()) {
            throw std::runtime_error("duplicated barcode sequences");
        } else {
            reference[curhash]=i;
        }
    }

    // Specifying the lengths that need to be tested.
    size_t nlens=0;
    for (auto v : valid_lengths) { nlens += v; }
    std::vector<int> all_lens;
    all_lens.reserve(nlens);
    for (size_t i=0; i<valid_lengths.size(); ++i) {
        if (valid_lengths[i]) { all_lens.push_back(i); }
    }

    return std::make_pair(reference, all_lens);
}
