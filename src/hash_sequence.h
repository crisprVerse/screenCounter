#ifndef SEQHASH_H
#define SEQHASH_H

#include <cstdint>
#include <vector>
#include <cstddef>

/* OVERVIEW:
 *
 * This file contains methods to convert a DNA sequence into a hash.
 * Each base can be stored in 2 bits, and we pack multiple bases 
 * into a 64-bit word. We then create a vector of these words
 * to enable hashing of barcodes that are longer than 32 bp.
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

#include <boost/functional/hash.hpp>

namespace std {

/* Defines a hashing function for the seqhash,
 * using Boost's hash_range function.
 */

template <>
struct hash<seqhash>{
    size_t operator()(const seqhash& k) const {
        return boost::hash_range(k.begin(), k.end());
    }
};

}

#endif
