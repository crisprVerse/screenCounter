#include "Rcpp.h"

#include "hash_sequence.h"
#include "build_dictionary.h"
#include "search_sequence.h"

#include <stdexcept>
#include <vector>
#include <sstream>
#include <map>

/* Combination guide parser. */

template<size_t nvariable>
struct result_store {
    std::vector<int> counts;

    // Only adding if all the indices agree with each other.
    void operator()(const combination<nvariable>& input) {
        if (nvariable==0) { return; }
        int ref=input.indices[0];
        for (size_t i=1; i<nvariable; ++i) {
            if (input.indices[i]!=ref) { return; }
        }
        ++(counts[ref]);
        return;
    }
};

template<size_t nvariable> 
struct se_fixed_info {
    se_fixed_info(Rcpp::List g, Rcpp::StringVector c, Rcpp::LogicalVector s, Rcpp::LogicalVector d) :
        info(g, c, s, d) 
    {
        if (info.nbarcodes.size()) {
            size_t ref=info.nbarcodes.front();
            for (size_t i=1; i<info.nbarcodes.size(); ++i) {
                if (ref!=info.nbarcodes[i]){ 
                    throw std::runtime_error("number of sequences differ between variable regions");
                }
            }
            output.counts.resize(ref);
        }
    }
    search_info<nvariable> info;
    result_store<nvariable> output;
};

template<size_t nvariable>
SEXP setup_barcodes_fixed(SEXP constants, SEXP guide_list, Rcpp::LogicalVector allowS, Rcpp::LogicalVector allowD) 
{
    se_fixed_info<nvariable> * ptr = new se_fixed_info<nvariable>(
        guide_list, constants, allowS, allowD
    );
    return Rcpp::XPtr<se_fixed_info<nvariable> >(ptr, true);
}

template<size_t nvariable>
SEXP count_barcodes_fixed(SEXP seqs, SEXP xptr) {
    Rcpp::XPtr<se_fixed_info<nvariable> > ptr(xptr);
    const auto& search_info=ptr->info;
    auto& output=ptr->output;

    Rcpp::StringVector Seqs(seqs);
    for (size_t i=0; i<Seqs.size(); ++i) {
        Rcpp::String s=Seqs[i];
        const char* sptr=s.get_cstring();
        const size_t len=Rf_length(s.get_sexp());
        search_sequence(sptr, len, search_info, output);
    }

    return R_NilValue;
}

template<size_t nvariable>
SEXP report_barcodes_fixed(SEXP xptr) {
    Rcpp::XPtr<se_fixed_info<nvariable> > ptr(xptr);
    const auto& out_store=ptr->output;
    return Rcpp::IntegerVector(out_store.counts.begin(), out_store.counts.end());
}

/*********************************************************
 * Realizations of template functions for 1 or 2 guides. *
 *********************************************************/

// [[Rcpp::export(rng=false)]]
SEXP setup_barcodes_fixed_solo(SEXP constants, SEXP guide_list, Rcpp::LogicalVector allowSub, Rcpp::LogicalVector allowDel) {
    return setup_barcodes_fixed<1>(constants, guide_list, allowSub, allowDel);
}

// [[Rcpp::export(rng=false)]]
SEXP count_barcodes_fixed_solo(SEXP seqs, SEXP xptr) {
    return count_barcodes_fixed<1>(seqs, xptr);
}

// [[Rcpp::export(rng=false)]]
SEXP report_barcodes_fixed_solo(SEXP xptr) {
    return report_barcodes_fixed<1>(xptr);
}

// [[Rcpp::export(rng=false)]]
SEXP setup_barcodes_fixed_dual(SEXP constants, SEXP guide_list, Rcpp::LogicalVector allowSub, Rcpp::LogicalVector allowDel) {
    return setup_barcodes_fixed<2>(constants, guide_list, allowSub, allowDel);
}

// [[Rcpp::export(rng=false)]]
SEXP count_barcodes_fixed_dual(SEXP seqs, SEXP xptr) {
    return count_barcodes_fixed<2>(seqs, xptr);
}

// [[Rcpp::export(rng=false)]]
SEXP report_barcodes_fixed_dual(SEXP xptr) {
    return report_barcodes_fixed<2>(xptr);
}

