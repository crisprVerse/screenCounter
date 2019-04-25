#ifndef HASH_SEQUENCE_H
#define HASH_SEQUENCE_H

#include <cstdint>
#include <vector>
#include <cstddef>

/* OVERVIEW:
 *
 * This file contains methods to convert a DNA sequence into a hash.
 * Each base can be stored in 2 bits, and we use a std::vector 
 * pack multiple bases into a 64-bit word. We use a string
 * in order to take advantage of the C++ standard's 
 * in-built support for string hashing in an unordered_map.
 *
 * We also implement classes to modify a pre-built hash for a 
 * given sequence, by performing substitutions or deletions.
 * This avoids the need to recompute the hash in its entirety
 * if we want to compute slight variations on a sequence.
 */

class seqhash {
public:    
    seqhash(size_t=0);
    typedef std::uint64_t word;
    bool operator==(const seqhash&) const;
    bool operator!=(const seqhash&) const;

    size_t size () const;
    std::vector<word>::iterator begin();
    std::vector<word>::iterator end();
    std::vector<word>::const_iterator begin() const;
    std::vector<word>::const_iterator end() const;

    word& operator[](size_t);
    const word& operator[](size_t) const;
    void pop_back();
private:    
    std::vector<word> words;
};

seqhash hash_sequence(const char*, const size_t);

void shift_sequence(seqhash&, const size_t, char);

struct rolling_substitution {
public:
    rolling_substitution(const seqhash&, size_t);
    bool advance();
    const seqhash& get() const;
private:
    seqhash hash;
    size_t len, word, pos;
    seqhash::word current, state, original;
};

struct rolling_deletion {
public:
    rolling_deletion(const seqhash&, size_t);
    bool advance();
    const seqhash& get() const;
private:
    seqhash hash;
    size_t word, pos;
    seqhash::word current, discarded;
};

#endif
