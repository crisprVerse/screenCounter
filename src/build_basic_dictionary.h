#ifndef SEQUENCE_DICTIONARY_H
#define SEQUENCE_DICTIONARY_H

#include "Rcpp.h"
#include "hash_sequence.h"
#include <unordered_map>

typedef std::unordered_map<seqhash, int> basic_dictionary;

basic_dictionary build_basic_dictionary(Rcpp::StringVector);

#endif
