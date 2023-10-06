#ifndef UTILS_H
#define UTILS_H

#include "Rcpp.h"
#include "kaori/kaori.hpp"
#include <vector>

kaori::BarcodePool format_pointers(const Rcpp::CharacterVector&);

kaori::SearchStrand to_strand(bool);

kaori::SearchStrand to_strand(int);

template<size_t V>
std::pair<Rcpp::IntegerMatrix, Rcpp::IntegerVector> count_combinations(const std::vector<std::array<int, V> >& available) {
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

    return std::make_pair(std::move(indices), std::move(freq));
}

#endif
