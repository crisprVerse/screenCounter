#include "Rcpp.h"

#include "kaori/kaori.hpp"
#include "byteme/byteme.hpp"

#include <algorithm>
#include <vector>

#include "utils.h"

template<size_t N, class Reader>
Rcpp::List count_dual_barcodes_single_end_(
    Reader& reader, 
    std::string constant, 
    const std::vector<kaori::BarcodePool>& pools, 
    int strand,
    int mismatches,
    bool use_first, 
    bool diagnostics,
    int nthreads)
{
    typename kaori::DualBarcodesSingleEnd<N>::Options options;
    options.strand = to_strand(strand);
    options.max_mismatches = mismatches;
    options.use_first = use_first;

    if (!diagnostics) {
        kaori::DualBarcodesSingleEnd<N> handler(constant.c_str(), constant.size(), pools, options);
        kaori::process_single_end_data(&reader, handler, nthreads);
        const auto& counts = handler.get_counts();
        return Rcpp::List::create(
            Rcpp::IntegerVector(counts.begin(), counts.end()),
            Rcpp::IntegerVector::create(handler.get_total())
        );

    } else {
        kaori::DualBarcodesSingleEndWithDiagnostics<N, 2> handler(constant.c_str(), constant.size(), pools, options);
        kaori::process_single_end_data(&reader, handler, nthreads);
        handler.sort();
        const auto& counts = handler.get_counts();
        return Rcpp::List::create(
            Rcpp::IntegerVector(counts.begin(), counts.end()),
            count_combinations(handler.get_combinations(), 3),
            Rcpp::IntegerVector::create(handler.get_total())
        );
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::List count_dual_barcodes_single_end(
    std::string path, 
    std::string constant, 
    Rcpp::List pools,
    int strand,
    int mismatches,
    bool use_first, 
    bool diagnostics,
    int nthreads)
{
    byteme::SomeFileReader reader(path.c_str());

    std::vector<kaori::BarcodePool> ptr_pools;
    ptr_pools.reserve(pools.size());
    for (size_t p = 0, end = pools.size(); p < end; ++p) {
        Rcpp::CharacterVector current(pools[p]);
        ptr_pools.push_back(format_pointers(current));
    }

    size_t len = constant.size();
    Rcpp::List output;
    if (len <= 32) {
        output = count_dual_barcodes_single_end_< 32>(reader, constant, ptr_pools, strand, mismatches, use_first, diagnostics, nthreads);
    } else if (len <= 64) {
        output = count_dual_barcodes_single_end_< 64>(reader, constant, ptr_pools, strand, mismatches, use_first, diagnostics, nthreads);
    } else if (len <= 128) {
        output = count_dual_barcodes_single_end_<128>(reader, constant, ptr_pools, strand, mismatches, use_first, diagnostics, nthreads);
    } else if (len <= 256) {
        output = count_dual_barcodes_single_end_<256>(reader, constant, ptr_pools, strand, mismatches, use_first, diagnostics, nthreads);
    } else {
        throw std::runtime_error("lacking compile-time support for constant regions longer than 256 bp");
    }

    return output;
}

