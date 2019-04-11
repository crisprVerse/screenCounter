#include "Rcpp.h"
#include <string>
#include <unordered_map>
#include <vector>

struct sequence_dictionary {
    std::unordered_map<std::u32string, std::pair<int, int> > mapping;
    size_t len=0;
};

sequence_dictionary build_dictionary(Rcpp::StringVector, bool);

sequence_dictionary build_deleted_dictionary(Rcpp::StringVector);

