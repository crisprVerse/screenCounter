#include "Rcpp.h"

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

    Rcpp::StringVector Seqs(seqs);
    for (size_t i=0; i<Seqs.size(); ++i) {
        Rcpp::String s=Seqs[i];
        const char* sptr=s.get_cstring();
        const size_t len=Rf_length(s.get_sexp());
        search_sequence<1>(sptr, len, search_info, output);
    }

    return R_NilValue;
}

// [[Rcpp::export(rng=false)]]
SEXP report_barcodes_single(SEXP xptr) {
    Rcpp::XPtr<se_single_info> ptr(xptr);
    const auto& out_store=ptr->output;
    return Rcpp::IntegerVector(out_store.counts.begin(), out_store.counts.end());
}

