#ifndef UTILS_H
#define UTILS_H

#include "Rcpp.h"

void reverse_complement (char*, size_t);

bool is_valid (char);

int is_valid (const char*, size_t);

size_t check_length (const Rcpp::StringVector&);

#endif
