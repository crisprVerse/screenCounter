#include "Rcpp.h"

extern "C" {

#include "Biostrings_interface.h"

}

#include "hash_sequence.h"
#include "build_dictionary.h"
#include "search_sequence.h"

#include <stdexcept>
#include <vector>
#include <sstream>
#include <map>

/* Combination guide parser. */

struct result_store {
    std::vector<int> counts;
    void operator()(const combination<1>& input) {
        ++(counts[input.indices[0]]);
        return;
    }
};

struct se_single_info {
    se_single_info(Rcpp::List g, Rcpp::StringVector c, Rcpp::LogicalVector s, Rcpp::LogicalVector d) :
        info(g, c, s, d)
    {
        size_t ref=info.nbarcodes.front();
        output.counts.resize(ref);
    }
    search_info<1> info;
    result_store output;
};

// [[Rcpp::export(rng=false)]]
SEXP setup_barcodes_single(SEXP constants, SEXP guide_list, Rcpp::LogicalVector allowS, Rcpp::LogicalVector allowD) 
{
    se_single_info * ptr = new se_single_info(guide_list, constants, allowS, allowD);
    return Rcpp::XPtr<se_single_info>(ptr, true);
}

// [[Rcpp::export(rng=false)]]
SEXP count_barcodes_single(SEXP seqs, SEXP xptr) {
    Rcpp::XPtr<se_single_info> ptr(xptr);
    const auto& search_info=ptr->info;
    auto& output=ptr->output;

    auto Seqs=hold_XStringSet(seqs);
    size_t nseqs=get_length_from_XStringSet_holder(&Seqs);
    std::vector<char> curseq;

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

        search_sequence<1>(curseq.data(), len, search_info, output);
    }

    return R_NilValue;
}

// [[Rcpp::export(rng=false)]]
SEXP report_barcodes_single(SEXP xptr) {
    Rcpp::XPtr<se_single_info> ptr(xptr);
    const auto& out_store=ptr->output;
    return Rcpp::IntegerVector(out_store.counts.begin(), out_store.counts.end());
}

