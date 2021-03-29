#include "Rcpp.h"
#include "hash_sequence.h"
#include "build_basic_dictionary.h"
#include "build_combined_dictionary.h"
#include <limits>

/* Test code to check that the hashing works correctly. */

Rcpp::NumericVector hash2double(const seqhash& input) {
    Rcpp::NumericVector output(input.size()*2);
    for (size_t i=0; i<input.size(); ++i) {
        auto current=input[i];

        // Get into a [0, 2^32) range to guarantee integer representation in a double.
        int SHIFT=32;
        auto upper_half = current >> SHIFT;
        auto lower_half = current & ((static_cast<seqhash::word>(1) << SHIFT) - 1);

        output[2*i] = static_cast<double>(lower_half);
        output[2*i+1] = static_cast<double>(upper_half);
    }
    return output;
}

// [[Rcpp::export(rng=false)]]
SEXP basic_hash(Rcpp::StringVector input) {
    std::string s=Rcpp::as<std::string>(input[0]);
    auto as_hash=hash_sequence(s.c_str(), s.size());
    return hash2double(as_hash);
}

// [[Rcpp::export(rng=false)]]
SEXP shift_hash(Rcpp::StringVector input, Rcpp::StringVector coming) {
    std::string s=Rcpp::as<std::string>(input[0]);
    auto as_hash=hash_sequence(s.c_str(), s.size());

    std::string next=Rcpp::as<std::string>(coming[0]);
    Rcpp::List output(next.size());
    for (size_t i=0; i<next.size(); ++i) {
        shift_sequence(as_hash, s.size(), next[i]);
        output[i]=hash2double(as_hash);
    }

    return output;
}

// [[Rcpp::export(rng=false)]]
SEXP export_basic_dictionary(Rcpp::StringVector guides) {
    auto mapping = build_basic_dictionary(guides);
    Rcpp::List keys(mapping.size());
    Rcpp::IntegerVector idx(mapping.size());

    size_t i=0;
    for (const auto& KV : mapping) {
        keys[i] = hash2double(KV.first);
        idx[i] = KV.second + 1; // get back to 1-based indices.
        ++i;
    }

    return Rcpp::List::create(keys, idx);
}

// [[Rcpp::export(rng=false)]]
SEXP export_combined_dictionary(Rcpp::StringVector constant, Rcpp::List variable, int n_sub, int n_insert, int n_del, int n_total) { 
    combined_dictionary cdict = build_combined_dictionary(constant, variable, n_sub, n_insert, n_del, n_total);

    Rcpp::List output(cdict.size());
    size_t c = 0;
    for (const auto& CV : cdict) {
        const auto& mapping = CV.second;
        Rcpp::List keys(mapping.size()), indices(mapping.size());
        Rcpp::IntegerVector priority(mapping.size());

        size_t i = 0;
        for (const auto& KV : mapping) {
            keys[i] = hash2double(KV.first);
            priority[i] = KV.second.first;

            const auto& v = KV.second.second;
            Rcpp::IntegerVector tmp (v.begin(), v.end());
            for (auto& j : tmp) { ++j; } // get back to 1-based indices.
            indices[i] = tmp;

            ++i;
        }

        output[c] = Rcpp::List::create(Rcpp::IntegerVector::create(CV.first), keys, priority, indices);
        ++c;
    }

    return output;
}

