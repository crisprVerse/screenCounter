#ifndef HASH_SEQUENCE_H
#define HASH_SEQUENCE_H
#include <string>

/* OVERVIEW:
 *
 * This file contains methods to convert a DNA sequence into a hash.
 * Each base can be stored in 2 bits, and we use a u32string 
 * pack multiple bases into a 32-bit word. We use a string
 * in order to take advantage of the C++ standard's 
 * in-built support for string hashing in an unordered_map.
 *
 * We also implement classes to modify a pre-built hash for a 
 * given sequence, by performing substitutions or deletions.
 * This avoids the need to recompute the hash in its entirety
 * if we want to compute slight variations on a sequence.
 */

std::u32string hash_sequence(const char*, const size_t);

void shift_sequence(std::u32string&, const size_t, char);

struct rolling_substitution {
public:
    rolling_substitution(const std::u32string&, size_t);
    bool advance();
    const std::u32string& get() const;
private:
    std::u32string hash;
    size_t len, word, pos;
    uint32_t current, state, original;
};

struct rolling_deletion {
public:
    rolling_deletion(const std::u32string&, size_t);
    bool advance();
    const std::u32string& get() const;
private:
    std::u32string hash;
    size_t word, pos;
    uint32_t current, discarded;
};

#endif
