#include "Rcpp.h"

#include "kaori/handlers/SingleBarcodeSingleEnd.hpp"
#include "kaori/process_data.hpp"
#include "byteme/SomeFileReader.hpp"

#include <algorithm>
#include <vector>

#include "utils.h"

template<size_t N, class Reader>
void count_single_barcodes_(Rcpp::IntegerVector& output, int& total_nreads, Reader& reader, std::string constant, int strand, const kaori::BarcodePool& options, int mismatches, bool use_first, int nthreads) {
    kaori::SingleBarcodeSingleEnd<N> handler(constant.c_str(), constant.size(), strand, options, mismatches);
    handler.set_first(use_first);
    kaori::process_single_end_data(&reader, handler, nthreads);

    const auto& counts = handler.get_counts();
    std::copy(counts.begin(), counts.end(), output.begin());
    total_nreads = handler.get_total();

    return;
}

//[[Rcpp::export(rng=false)]]
Rcpp::List count_single_barcodes(std::string path, std::string constant, int strand, Rcpp::CharacterVector options, int mismatches, bool use_first, int nthreads) {
    byteme::SomeFileReader reader(path.c_str());
    auto ptrs = format_pointers(options);

    Rcpp::IntegerVector output_counts(options.size());
    int output_totals;

    // Support up to 256 bp constant regions.
    if (constant.size() <= 32) {
        count_single_barcodes_<32>(output_counts, output_totals, reader, constant, strand, ptrs, mismatches, use_first, nthreads);
    } else if (constant.size() <= 64) {
        count_single_barcodes_<64>(output_counts, output_totals, reader, constant, strand, ptrs, mismatches, use_first, nthreads);
    } else if (constant.size() <= 128) {
        count_single_barcodes_<32>(output_counts, output_totals, reader, constant, strand, ptrs, mismatches, use_first, nthreads);
    } else if (constant.size() <= 256) {
        count_single_barcodes_<64>(output_counts, output_totals, reader, constant, strand, ptrs, mismatches, use_first, nthreads);
    } else {
        throw std::runtime_error("lacking compile-time support for constant regions longer than 256 bp");
    }

    return Rcpp::List::create(output_counts, output_totals);
}
