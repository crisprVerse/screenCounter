#include "Rcpp.h"

#include "kaori/handlers/CombinatorialBarcodes.hpp"
#include "kaori/process_data.hpp"
#include "byteme/SomeFileReader.hpp"

#include <algorithm>
#include <vector>

#include "utils.h"

template<size_t N, size_t V, class Reader>
Rcpp::List count_combo_barcodes_(Reader& reader, 
    std::string constant, 
    int strand, 
    const std::vector<std::vector<const char*> >& options,
    int mismatches, 
    bool use_first, 
    int nthreads) 
{
    kaori::CombinatorialBarcodes<N, V> handler(constant.c_str(), constant.size(), strand, options, mismatches);
    handler.set_first(use_first);
    kaori::process_single_end_data(&reader, handler, nthreads);

    handler.sort();
    const auto& available = handler.get_combinations();

    std::vector<size_t> boundaries;
    if (available.size() > 0) {
        boundaries.push_back(0);
        for (size_t i = 1; i < available.size(); ++i) {
            if (available[i] != available[i-1]) {
                boundaries.push_back(i);
            }
        }
    }

    auto count = boundaries.size();

    Rcpp::IntegerMatrix indices(V, count);
    auto it = indices.begin();
    for (auto b : boundaries) {
        const auto& current = available[b];
        std::copy(current.begin(), current.end(), it);
        it += V;
    }

    Rcpp::IntegerVector freq(count);
    if (count) {
        for (size_t b = 1; b < count; ++b) {
            freq[b - 1] = boundaries[b] - boundaries[b - 1];
        }
        freq[count - 1] = available.size() - boundaries[count - 1];
    }

    return Rcpp::List::create(indices, freq, Rcpp::IntegerVector::create(handler.get_total()));
}

//[[Rcpp::export(rng=false)]]
Rcpp::List count_combo_barcodes(std::string path, std::string constant, int strand, Rcpp::List options, int mismatches, bool use_first, int nthreads) {
    byteme::SomeFileReader reader(path.c_str());

    std::vector<std::vector<const char*> > opts;
    std::vector<Rcpp::CharacterVector> converted;
    for (size_t o = 0; o < options.size(); ++o) {
        converted.push_back(Rcpp::CharacterVector(options[o]));
        opts.push_back(format_pointers(converted.back()));
    }

    // Support up to 256 bp constant regions.
    Rcpp::List output;
    if (constant.size() <= 32) {
        output = count_combo_barcodes_<128, 2>(reader, constant, strand, opts, mismatches, use_first, nthreads);
    } else if (constant.size() <= 64) {
        output = count_combo_barcodes_<256, 2>(reader, constant, strand, opts, mismatches, use_first, nthreads);
    } else if (constant.size() <= 128) {
        output = count_combo_barcodes_<512, 2>(reader, constant, strand, opts, mismatches, use_first, nthreads);
    } else if (constant.size() <= 256) {
        output = count_combo_barcodes_<1024, 2>(reader, constant, strand, opts, mismatches, use_first, nthreads);
    } else {
        throw std::runtime_error("lacking compile-time support for constant regions longer than 256 bp");
    }

    return output;
}
