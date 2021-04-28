#include "Rcpp.h"

extern "C" {

#include "Biostrings_interface.h"

}

#include "hash_sequence.h"
#include "build_combined_dictionary.h"
#include "simple_search.h"
#include "utils.h"

#include <vector>

struct pe_dual_info {
    pe_dual_info(const Rcpp::StringVector& constants, const Rcpp::StringVector& variables, int n_sub, int n_insert, int n_del, int n_total) {
        std::vector<Rcpp::StringVector> varlist(1);
        varlist[0] = variables;
        mapping = build_combined_dictionary(constants, varlist, n_sub, n_insert, n_del, n_total);
        return;
    }
    combined_dictionary mapping;
};

// [[Rcpp::export(rng=false)]]
SEXP setup_barcodes_dual(Rcpp::StringVector constants, Rcpp::StringVector variables, int n_sub, int n_insert, int n_del, int n_total) {
    pe_dual_info * ptr = new pe_dual_info(constants, variables, n_sub, n_insert, n_del, n_total);
    return Rcpp::XPtr<pe_dual_info>(ptr, true);
}

// [[Rcpp::export(rng=false)]]
Rcpp::List count_barcodes_dual(SEXP seqs1, SEXP seqs2, SEXP xptr1, SEXP xptr2, bool forward1, bool forward2, bool randomized) {
    Rcpp::XPtr<pe_dual_info> ptr1(xptr1);
    const auto& edict1 = ptr1->mapping;
    Rcpp::XPtr<pe_dual_info> ptr2(xptr2);
    const auto& edict2 = ptr2->mapping;

    auto Seqs1 = hold_XStringSet(seqs1);
    size_t nseqs = get_length_from_XStringSet_holder(&Seqs1);
    auto Seqs2 = hold_XStringSet(seqs2);
    if (nseqs != get_length_from_XStringSet_holder(&Seqs2)) {
        throw std::runtime_error("paired FASTQ files should have the same number of sequences");
    }
    std::vector<char> curseq;

    Rcpp::IntegerVector out1(nseqs), out2(nseqs), diagnostics(2);
    int& nprovided = diagnostics[0];
    int& nother = diagnostics[1];
    constexpr int UNMAPPED = -1;

    for (size_t i=0; i<nseqs; ++i) {
        int id1_provided=UNMAPPED, id2_provided=UNMAPPED, id1_other=UNMAPPED, id2_other=UNMAPPED;
        int provided_nedits=0, other_nedits=0;

        for (int j=0; j<2; ++j) {
            auto sptr = (j==0 ? &Seqs1 : &Seqs2);

            // Search parameters for the provided orientation, i.e., read1 = barcode1, read2 = barcode2.
            auto& edict_provided = (j==0 ? edict1 : edict2);
            int& id_provided = (j==0 ? id1_provided : id2_provided);
            const bool& forward_provided = (j==0 ? forward1 : forward2);

            // Search parameters for the other orientation, i.e., read1 = barcode2, read2 = barcode1.
            auto& edict_other = (j==0 ? edict2 : edict1);
            int& id_other = (j==0 ? id2_other : id1_other);
            const bool& forward_other = (j==0 ? forward2 : forward1);

            auto current = get_elt_from_XStringSet_holder(sptr, i);
            const char* ptr = current.ptr;
            const size_t len = current.length;

            if (len > curseq.size()) {
                curseq.resize(len+1);
            }
            for (size_t i=0; i<len; ++i) {
                curseq[i]=DNAdecode(ptr[i]);
            }

            /* Searching one or both strands. Note that the different strands
             * are never competing for the same identifier here, so there's no
             * need to carry over the 'best_nedits' and 'chosen' from one
             * strand to another. (In fact, doing so would be incorrect).
             */
            if (forward_provided) {
                int best_nedits = len; // any large value will do here.
                simple_search(edict_provided, curseq.data(), len, best_nedits, id_provided);
                provided_nedits += best_nedits;
            }

            if (randomized && forward_other) {
                int best_nedits = len; 
                simple_search(edict_other, curseq.data(), len, best_nedits, id_other);
                other_nedits += best_nedits;
            }

            if (!forward_provided || (randomized && !forward_other)) {
                reverse_complement(curseq.data(), len);

                if (!forward_provided) {
                    int best_nedits = len; 
                    simple_search(edict_provided, curseq.data(), len, best_nedits, id_provided);
                    provided_nedits += best_nedits;
                }
                if (randomized && !forward_other) {
                    int best_nedits = len; 
                    simple_search(edict_other, curseq.data(), len, best_nedits, id_other);
                    other_nedits += best_nedits;
                }
            }
        }

        if (!randomized) {
            out1[i] = id1_provided;
            out2[i] = id2_provided;
        } else {
            const bool provided_ok = id1_provided != UNMAPPED && id2_provided != UNMAPPED;
            const bool other_ok = id1_other != UNMAPPED && id2_other != UNMAPPED;

            if (!other_ok) {
                out1[i] = id1_provided;
                out2[i] = id2_provided;
                nprovided += provided_ok;
            } else if (!provided_ok) {
                out1[i] = id1_other;
                out2[i] = id2_other;
                nother += other_ok;
            } else {
                if (provided_nedits < other_nedits) {
                    out1[i] = id1_provided;
                    out2[i] = id2_provided;
                    ++nprovided;
                } else if (provided_nedits > other_nedits) {
                    out1[i] = id1_other;
                    out2[i] = id2_other;
                    ++nother;
                } else if (provided_nedits == other_nedits && id1_provided == id1_other && id2_provided == id2_other) {
                    out1[i] = id1_provided;
                    out2[i] = id2_provided;
                    ++nprovided;
                } else {
                    // Ties are considered unmapped.
                    out1[i] = -1;
                    out2[i] = -1;
                }
            }
        }

        // Getting back to 1-based indexing.
        ++out1[i];
        ++out2[i];
    }

    return Rcpp::List::create(out1, out2, diagnostics);
}
