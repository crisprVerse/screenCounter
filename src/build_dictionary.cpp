#include "build_dictionary.h"
#include "hash_sequence.h"

/* Builds a hash index, returning it alongside a vector of lengths.
 * The latter is necessary to accommodate variable length barcodes.
 * It is okay to assume that there are no N's here, as any N's 
 * should have been eliminated by barcode construction at the R level.
 */

sequence_dictionary build_dictionary (Rcpp::StringVector guides, bool allowSub) {
    sequence_dictionary reference;
    auto& mapping=reference.mapping;
    size_t& reflen=reference.len;

    for (size_t i=0; i<guides.size(); ++i) {
        Rcpp::String g=guides[i];
        const size_t len=Rf_length(g.get_sexp());

        // Checking that the length exists.
        if (reflen==0) {
            reflen=len;
        } else if (reflen!=len) {
            throw std::runtime_error("barcodes are of variable length");
        }

        // Adding the true sequence, which is a high priority (0).
        const char* ptr=g.get_cstring();
        auto curhash=hash_sequence(ptr, len);

        auto it=mapping.find(curhash);
        if (it!=mapping.end() && (it->second).second==0) {
            throw std::runtime_error("duplicated barcode sequences");
        } else {
            mapping[curhash]=std::make_pair(i, 0);
        }

        // Adding a sequence with one substitution, which is lower priority (1).
        if (allowSub) {
            rolling_substitution Subber(curhash, len);
            do {
                const auto& subhash=Subber.get();
                auto it=mapping.find(subhash);
                if (it!=mapping.end()) {
                    if ((it->second).second==1 && (it->second).first!=i) {
                        (it->second).first=-1; // indicating that it clashes with another mismatch.
                    }
                } else {
                    mapping[subhash]=std::make_pair(i, 1);
                }
            } while (Subber.advance());
        }
    }

    return reference;
}

/* Builds a hash index where one base has been deleted.
 * There's no need to worry about clashes with the full length
 * sequence because that's of a different length anyway and will never clash.
 */

sequence_dictionary build_deleted_dictionary (Rcpp::StringVector guides) {
    sequence_dictionary reference;
    auto& mapping=reference.mapping;
    size_t& reflen=reference.len;

    for (size_t i=0; i<guides.size(); ++i) {
        Rcpp::String g=guides[i];
        const size_t len=Rf_length(g.get_sexp());

        // Checking that the length exists.
        if (reflen==0) {
            reflen=len;
        } else if (reflen!=len) {
            throw std::runtime_error("barcodes are of variable length");
        }

        // Rolling through deletions to the current sequence.
        const char* ptr=g.get_cstring();
        auto curhash=hash_sequence(ptr, len);
        rolling_deletion Del(curhash, len);

        do {
            const auto& subhash=Del.get();
            auto it=mapping.find(subhash);
            if (it!=mapping.end()) {
                if ((it->second).second==1 && (it->second).first!=i) {
                    (it->second).first=-1; // indicating that it clashes with another deletion.
                }
            } else {
                mapping[subhash]=std::make_pair(i, 1);
            }
        } while(Del.advance());
    }

    --reflen;
    return reference;
}
