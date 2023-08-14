#include "Rcpp.h"

#include "kaori/kaori.hpp"
#include "utils.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List match_barcodes(Rcpp::CharacterVector sequences, Rcpp::CharacterVector choices, int substitutions, bool reverse) {
    typename kaori::SimpleBarcodeSearch::Options opt;
    opt.max_mismatches = substitutions;
    opt.reverse = reverse;

    auto pool = format_pointers(choices);
    kaori::SimpleBarcodeSearch searcher(pool, opt);
    auto state = searcher.initialize();

    Rcpp::IntegerVector output_id(sequences.size()), output_mm(sequences.size());
    auto oiIt = output_id.begin();
    auto omIt = output_mm.begin();

    auto x = format_pointers(sequences);
    for (auto p : x.pool) {
        searcher.search(p, state);

        if (state.index >= 0) {
            *oiIt = state.index + 1;
            *omIt = state.mismatches;
        } else {
            *oiIt = NA_INTEGER;
            *omIt = NA_INTEGER;
        }

        ++oiIt;
        ++omIt;
    }

    return Rcpp::List::create(output_id, output_mm);
}
