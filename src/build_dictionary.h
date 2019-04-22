#ifndef SEQUENCE_DICTIONARY_H
#define SEQUENCE_DICTIONARY_H

#include "Rcpp.h"
#include "hash_sequence.h"

#include <unordered_map>
#include <boost/functional/hash.hpp>

/* Defines a hashing function for the seqhash,
 * using Boost's hash_range function.
 */

namespace std {

template <>
struct hash<seqhash>{
    size_t operator()(const seqhash& k) const {
        return boost::hash_range(k.begin(), k.end());
    }
};

}

/* Defining the sequence dictionary class */

struct sequence_dictionary {
    std::unordered_map<seqhash, std::pair<int, int> > mapping;
    size_t len=0;
};

sequence_dictionary build_dictionary(Rcpp::StringVector, bool);

sequence_dictionary build_deleted_dictionary(Rcpp::StringVector);

#endif
