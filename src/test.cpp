#include "Rcpp.h"
#include "hash_sequence.h"
#include "build_dictionary.h"
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
SEXP substitute_hash(Rcpp::StringVector input) {
    std::string s=Rcpp::as<std::string>(input[0]);
    auto as_hash=hash_sequence(s.c_str(), s.size());
    rolling_substitution Subber(as_hash, s.size());

    Rcpp::List output(s.size()*3);
    size_t i=0;
    do {
        if (i==output.size()) {
            throw std::runtime_error("more substituted strings than expected");
        }
        output[i]=hash2double(Subber.get());
        ++i;
    } while (Subber.advance());

    return output;
}

// [[Rcpp::export(rng=false)]]
SEXP delete_hash(Rcpp::StringVector input) {
    std::string s=Rcpp::as<std::string>(input[0]);
    auto as_hash=hash_sequence(s.c_str(), s.size());
    rolling_deletion Delly(as_hash, s.size());

    Rcpp::List output(s.size());
    size_t i=0;
    do {
        if (i==output.size()) {
            throw std::runtime_error("more deleted strings than expected");
        }
        output[i]=hash2double(Delly.get());
        ++i;
    } while (Delly.advance());

    return output;
}

// [[Rcpp::export(rng=false)]]
SEXP build_dict(Rcpp::StringVector guides, bool allowS, bool allowD) {
    sequence_dictionary dict;
    if (allowD) {
        dict=build_deleted_dictionary(guides);
    } else {
        dict=build_dictionary(guides, allowS);
    }

    const auto& mapping=dict.mapping;
    Rcpp::List keys(mapping.size());
    Rcpp::IntegerVector idx(mapping.size()), priority(mapping.size());

    size_t i=0;
    for (const auto& KV : mapping) {
        keys[i]=hash2double(KV.first);
        idx[i]=KV.second.first+1; // get back to 1-based indices.
        priority[i]=KV.second.second;
        ++i;
    }

    return Rcpp::List::create(
        Rcpp::List::create(
            keys,
            Rcpp::List::create(idx, priority)
        ),
        Rcpp::IntegerVector::create(dict.len)
    );
}

