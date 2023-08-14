#include "Rcpp.h"

#include "kaori/kaori.hpp"
#include "byteme/byteme.hpp"

#include <algorithm>
#include <vector>

#include "utils.h"

template<size_t N, class Reader>
void count_random_barcodes_(Rcpp::List& output, int& total_nreads, Reader& reader, const std::string& constant, int strand, int mismatches, bool use_first, int nthreads) {
    typename kaori::RandomBarcodeSingleEnd<N>::Options options;
    options.strand = to_strand(strand);
    options.max_mismatches = mismatches;
    options.use_first = use_first;

    kaori::RandomBarcodeSingleEnd<N> handler(constant.c_str(), constant.size(), options);
    kaori::process_single_end_data(&reader, handler, nthreads);

    const auto& counts = handler.get_counts();
    Rcpp::CharacterVector sequences(counts.size());
    Rcpp::IntegerVector frequencies(counts.size());

    auto sIt = sequences.begin();
    auto fIt = frequencies.begin();
    for (const auto& pair : counts) {
        *sIt = pair.first;
        *fIt = pair.second;
        ++sIt;
        ++fIt;
    }

    output[0] = sequences;
    output[1] = frequencies;
    total_nreads = handler.get_total();
    return;
}

//[[Rcpp::export(rng=false)]]
Rcpp::List count_random_barcodes(std::string path, std::string constant, int strand, int mismatches, bool use_first, int nthreads) {
    byteme::SomeFileReader reader(path.c_str());
    Rcpp::List output_freq(2);
    int output_totals;

    // Support up to 256 bp constant regions.
    if (constant.size() <= 32) {
        count_random_barcodes_<32>(output_freq, output_totals, reader, constant, strand, mismatches, use_first, nthreads);
    } else if (constant.size() <= 64) {
        count_random_barcodes_<64>(output_freq, output_totals, reader, constant, strand, mismatches, use_first, nthreads);
    } else if (constant.size() <= 128) {
        count_random_barcodes_<128>(output_freq, output_totals, reader, constant, strand, mismatches, use_first, nthreads);
    } else if (constant.size() <= 256) {
        count_random_barcodes_<256>(output_freq, output_totals, reader, constant, strand, mismatches, use_first, nthreads);
    } else {
        throw std::runtime_error("lacking compile-time support for constant regions longer than 256 bp");
    }

    return Rcpp::List::create(output_freq, output_totals);
}
