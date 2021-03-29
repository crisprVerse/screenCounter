#ifndef BUILD_COMBINED_DICTIONARY_H
#define BUILD_COMBINED_DICTIONARY_H

#include "Rcpp.h"
#include "hash_sequence.h"
#include <unordered_map>
#include <vector>

typedef std::unordered_map<seqhash, std::pair<int, std::vector<int> > > combined_dictionary_part;

typedef std::unordered_map<int, combined_dictionary_part> combined_dictionary;

combined_dictionary build_combined_dictionary(const Rcpp::StringVector&, const Rcpp::List&, int, int, int, int);

combined_dictionary build_combined_dictionary(const Rcpp::StringVector&, const std::vector<Rcpp::StringVector>&, int, int, int, int);

#endif
