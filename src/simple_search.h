#ifndef SIMPLE_SEARCH_H
#define SIMPLE_SEARCH_H

#include "build_combined_dictionary.h"

void simple_search(const combined_dictionary&, const char*, size_t, int&, int&);

void simple_search(const combined_dictionary&, const char*, size_t, int&, std::vector<int>&);

#endif

