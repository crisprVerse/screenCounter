#include "Rcpp.h"

extern "C" {

#include "Biostrings_interface.h"

}

#include "hash_sequence.h"
#include "build_combined_dictionary.h"
#include "simple_search.h"
#include "utils.h"

#include <vector>

struct se_single_info {
    se_single_info(const Rcpp::StringVector& constants, const Rcpp::StringVector& variables, int n_sub, int n_insert, int n_del, int n_total) : nchoices(variables.size()) {
        std::vector<Rcpp::StringVector> varlist(1);
        varlist[0] = variables;
        mapping = build_combined_dictionary(constants, varlist, n_sub, n_insert, n_del, n_total);
        return;
    }

    int nchoices;
    combined_dictionary mapping;
};

// [[Rcpp::export(rng=false)]]
SEXP setup_barcodes_single(Rcpp::StringVector constants, Rcpp::StringVector variables, int n_sub, int n_insert, int n_del, int n_total) {
    se_single_info * ptr = new se_single_info(constants, variables, n_sub, n_insert, n_del, n_total);
    return Rcpp::XPtr<se_single_info>(ptr, true);
}

/* This template function assumes that the STORE class implements:
 *
 * - A constructor accepting the number of barcodes and the number of sequences.
 * - A () method that fulfills the requirements of search_sequence<1>.
 * - An advance() method that is run in the loop body.
 * - A yield() method to return an IntegerVector back to R.
 */

template<class STORE>
Rcpp::IntegerVector search_barcodes_single(SEXP seqs, SEXP xptr, bool use_forward, bool use_reverse) {
    Rcpp::XPtr<se_single_info> ptr(xptr);
    const auto& edict = ptr->mapping;

    auto Seqs=hold_XStringSet(seqs);
    size_t nseqs=get_length_from_XStringSet_holder(&Seqs);
    std::vector<char> curseq;

    STORE output(ptr->nchoices, nseqs);

    for (size_t i=0; i<nseqs; ++i) {
        auto current=get_elt_from_XStringSet_holder(&Seqs, i);
        const char* ptr=current.ptr;
        const size_t len=current.length;

        if (len > curseq.size()) {
            curseq.resize(len+1);
        }
        for (size_t i=0; i<len; ++i) {
            curseq[i]=DNAdecode(ptr[i]);
        }

        int best_nedits = len; // any large value will do here.
        int chosen = -1;

        // Searching one or both strands.
        if (use_forward) {
            simple_search(edict, curseq.data(), len, best_nedits, chosen);
        }

        if (best_nedits > 0 && use_reverse) {
            reverse_complement(curseq.data(), len);
            simple_search(edict, curseq.data(), len, best_nedits, chosen);
        }

        if (chosen >= 0) {
            output(chosen);
        }
        output.advance();
    }

    return output.yield();
}

/*** Counting the instances of each barcode ***/

struct result_count_store {
    result_count_store(int nbarcodes, int nsequences) : counts(nbarcodes) {}
    void operator()(int input) {
        ++(counts[input]);
        return;
    }
    void advance() {}
    Rcpp::IntegerVector yield() const {
        return counts;
    }
    Rcpp::IntegerVector counts;
};

// [[Rcpp::export(rng=false)]]
Rcpp::IntegerVector count_barcodes_single(SEXP seqs, SEXP xptr, bool use_forward, bool use_reverse) {
    return search_barcodes_single<result_count_store>(seqs, xptr, use_forward, use_reverse);
} 

/*** Listing the barcode for each sequence ***/

struct result_identity_store {
    result_identity_store(int nbarcodes, int nsequences) : ids(nsequences) {}
    void operator()(int input) {
        ids[counter]=input+1;
        return;
    }
    void advance() { 
        ++counter;
        return;
    }
    Rcpp::IntegerVector yield() const {
        return ids;
    }
    Rcpp::IntegerVector ids;
    int counter=0;
};

// [[Rcpp::export(rng=false)]]
Rcpp::IntegerVector identify_barcodes_single(SEXP seqs, SEXP xptr, bool use_forward, bool use_reverse) {
    return search_barcodes_single<result_identity_store>(seqs, xptr, use_forward, use_reverse);
} 
