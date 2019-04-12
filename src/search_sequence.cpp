#include "search_sequence.h"

/* Utility functions to check that a string or character is valid,
 * i.e., contains only ACTG (or lower case) bases.
 */

bool is_valid(char base) {
    switch (base) {
        case 'A': case 'a': case 'C': case 'c': case 'G': case 'g': case 'T': case 't':
            return true;
    };
    return false;
}

int is_valid(const char* ptr, size_t n) {
    int valid=0;
    for (size_t i=0; i<n; ++i) {
        valid+=is_valid(ptr[i]);
    }
    return valid;
}

/* Class that combines shift_sequence and is_valid to iterate across a
 * string and produce the hash at each position (or die trying).
 */

hash_scanner::hash_scanner(const char* p, size_t n) : ptr(p), len(n), 
    hashed(hash_sequence(p, n)), nvalid(is_valid(p, n)) {}

void hash_scanner::advance() {
    nvalid-=is_valid(*ptr);
    const char next=*(ptr+len);
    ++ptr;

    shift_sequence(hashed, len, next);
    nvalid+=is_valid(next);
    return;
}

const std::u32string& hash_scanner::hash() const { return hashed; }

bool hash_scanner::valid() const { return nvalid==len; }
