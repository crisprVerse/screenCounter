#include "Rcpp.h"
#include <string>
#include <unordered_map>
#include <vector>

typedef std::unordered_map<std::u32string, int> seqhash;

std::pair<seqhash, std::vector<int> > build_hash(Rcpp::StringVector);
