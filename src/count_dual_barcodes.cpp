#include "Rcpp.h"

#include "kaori/kaori.hpp"
#include "byteme/byteme.hpp"

#include <algorithm>
#include <vector>

#include "utils.h"

template<size_t N, class Reader>
Rcpp::List count_dual_barcodes_(
    Reader& reader1, 
    std::string constant1, 
    bool reverse1,
    const kaori::BarcodePool& pool1, 
    int mismatches1,

    Reader& reader2, 
    std::string constant2, 
    bool reverse2,
    const kaori::BarcodePool& pool2, 
    int mismatches2,

    bool randomized,
    bool use_first, 
    int nthreads) 
{
    typename kaori::DualBarcodes<N>::Options options;
    options.strand1 = to_strand(reverse1);
    options.max_mismatches1 = mismatches1;
    options.strand2 = to_strand(reverse2);
    options.max_mismatches2 = mismatches2;
    options.random = randomized;
    options.use_first = use_first;

    kaori::DualBarcodes<N> handler(
        constant1.c_str(), constant1.size(), pool1,
        constant2.c_str(), constant2.size(), pool2, 
        options
    );
    kaori::process_paired_end_data(&reader1, &reader2, handler, nthreads);

    const auto& counts = handler.get_counts();
    return Rcpp::List::create(
        Rcpp::IntegerVector(counts.begin(), counts.end()),
        Rcpp::IntegerVector::create(handler.get_total())
    );
}

template<size_t N, class Reader>
Rcpp::List count_dual_barcodes_diagnostics_(
    Reader& reader1, 
    std::string constant1, 
    bool reverse1,
    const kaori::BarcodePool& pool1, 
    int mismatches1,

    Reader& reader2, 
    std::string constant2, 
    bool reverse2,
    const kaori::BarcodePool& pool2, 
    int mismatches2,

    bool randomized,
    bool use_first, 
    int nthreads) 
{
    typename kaori::DualBarcodes<N>::Options options;
    options.strand1 = to_strand(reverse1);
    options.max_mismatches1 = mismatches1;
    options.strand2 = to_strand(reverse2);
    options.max_mismatches2 = mismatches2;
    options.random = randomized;
    options.use_first = use_first;

    kaori::DualBarcodesWithDiagnostics<N> handler(
        constant1.c_str(), constant1.size(), pool1, 
        constant2.c_str(), constant2.size(), pool2,
        options
    );
    kaori::process_paired_end_data(&reader1, &reader2, handler, nthreads);
    handler.sort();

    const auto& counts = handler.get_counts();
    return Rcpp::List::create(
        Rcpp::IntegerVector(counts.begin(), counts.end()),
        count_combinations(handler.get_combinations(), 3),
        Rcpp::IntegerVector::create(handler.get_total()),
        Rcpp::IntegerVector::create(handler.get_barcode1_only()),
        Rcpp::IntegerVector::create(handler.get_barcode2_only())
    );
}

//[[Rcpp::export(rng=false)]]
Rcpp::List count_dual_barcodes(
    std::string path1, 
    std::string constant1, 
    bool reverse1,
    int mismatches1,
    Rcpp::CharacterVector pool1,
        
    std::string path2, 
    std::string constant2, 
    bool reverse2,
    int mismatches2,
    Rcpp::CharacterVector pool2,

    bool randomized,
    bool use_first, 
    bool diagnostics,
    int nthreads)
{
    byteme::SomeFileReader reader1(path1.c_str());
    auto ptrs1 = format_pointers(pool1);

    byteme::SomeFileReader reader2(path2.c_str());
    auto ptrs2 = format_pointers(pool2);

    size_t len = std::max(constant1.size(), constant2.size());
    Rcpp::List output;

    if (!diagnostics) {
        // Support up to 256 bp constant regions.
        if (len <= 32) {
            output = count_dual_barcodes_<32>(reader1, constant1, reverse1, ptrs1, mismatches1, reader2, constant2, reverse2, ptrs2, mismatches2, randomized, use_first, nthreads);
        } else if (len <= 64) {
            output = count_dual_barcodes_<64>(reader1, constant1, reverse1, ptrs1, mismatches1, reader2, constant2, reverse2, ptrs2, mismatches2, randomized, use_first, nthreads);
        } else if (len <= 128) {
            output = count_dual_barcodes_<128>(reader1, constant1, reverse1, ptrs1, mismatches1, reader2, constant2, reverse2, ptrs2, mismatches2, randomized, use_first, nthreads);
        } else if (len <= 256) {
            output = count_dual_barcodes_<256>(reader1, constant1, reverse1, ptrs1, mismatches1, reader2, constant2, reverse2, ptrs2, mismatches2, randomized, use_first, nthreads);
        } else {
            throw std::runtime_error("lacking compile-time support for constant regions longer than 256 bp");
        }
    } else {
        if (len <= 32) {
            output = count_dual_barcodes_diagnostics_<32>(reader1, constant1, reverse1, ptrs1, mismatches1, reader2, constant2, reverse2, ptrs2, mismatches2, randomized, use_first, nthreads);
        } else if (len <= 64) {
            output = count_dual_barcodes_diagnostics_<64>(reader1, constant1, reverse1, ptrs1, mismatches1, reader2, constant2, reverse2, ptrs2, mismatches2, randomized, use_first, nthreads);
        } else if (len <= 128) {
            output = count_dual_barcodes_diagnostics_<128>(reader1, constant1, reverse1, ptrs1, mismatches1, reader2, constant2, reverse2, ptrs2, mismatches2, randomized, use_first, nthreads);
        } else if (len <= 256) {
            output = count_dual_barcodes_diagnostics_<256>(reader1, constant1, reverse1, ptrs1, mismatches1, reader2, constant2, reverse2, ptrs2, mismatches2, randomized, use_first, nthreads);
        } else {
            throw std::runtime_error("lacking compile-time support for constant regions longer than 256 bp");
        }
    }

    return output;
}

