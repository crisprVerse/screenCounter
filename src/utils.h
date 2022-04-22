#ifndef UTILS_H
#define UTILS_H

#include "Rcpp.h"
#include <vector>

std::vector<const char*> format_pointers(const Rcpp::CharacterVector&);

template<size_t V>
Rcpp::List count_combinations(const std::vector<std::array<int, V> >& available, int extras) {
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

    Rcpp::List output(extras + 2);
    output[0] = indices;
    output[1] = freq;
    return output;
}

#endif
