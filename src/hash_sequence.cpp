#include "hash_sequence.h"
#include <stdexcept>
#include <limits>

/*****************************
 **** seqhash definitions ****
 *****************************/

seqhash::seqhash(size_t n) : words(n) {}

size_t seqhash::size () const { return words.size(); }

std::vector<seqhash::word>::const_iterator seqhash::begin() const { return words.begin(); }

std::vector<seqhash::word>::iterator seqhash::begin() { return words.begin(); }

std::vector<seqhash::word>::const_iterator seqhash::end() const { return words.end(); }

std::vector<seqhash::word>::iterator seqhash::end() { return words.end(); }

void seqhash::pop_back() {
    words.pop_back();
    return;
}

bool seqhash::operator==(const seqhash& other) const {
    return words==other.words;
}

bool seqhash::operator!=(const seqhash& other) const {
    return words!=other.words;
}

seqhash::word& seqhash::operator[](size_t i) {
    return words[i];
}

const seqhash::word& seqhash::operator[](size_t i) const {
    return words[i];
}

#ifdef DEBUG
#include <iostream>
#include <bitset>
void disgorge(const seqhash& in) {
    typedef std::bitset<std::numeric_limits<seqhash::word>::digits> bitter;

    for (size_t i=0; i<in.size(); ++i) {
        std::cout << bitter(in[i]) << ", ";
    }
    std::cout << std::endl;
    return;
}

#endif

/*********************************
 **** basic hashing functions ****
 *********************************/

seqhash::word bitify (char base) {
    switch (base) {
        case 'A': case 'a':
            return 0;
        case 'C': case 'c':
            return 1;
        case 'G': case 'g':
            return 2;
        case 'T': case 't':
            return 3;
    };
    return 0; // N's and other weird crap goes here.
}

/* Converts a C-style string into a compressed string, 
 * for straightforward hashing. In each char, the two
 * least significant bits hold the 5'-most base.
 */

constexpr int BITS_PER_BASE=2;
constexpr int WIDTH_IN_BASES=std::numeric_limits<seqhash::word>::digits/BITS_PER_BASE;

seqhash hash_sequence(const char* p, const size_t n) { 
    int effective_length=(n + WIDTH_IN_BASES - 1)/WIDTH_IN_BASES;
    seqhash output(effective_length);

    size_t i=0, j=0;
    while (i < n) {
        size_t limit=std::min(n, i+WIDTH_IN_BASES);

        seqhash::word& current=output[j++];
        int shift=0;
        for (; i<limit; ++i) {
            current+=(bitify(p[i]) << shift);
            shift+=BITS_PER_BASE;
        }
    }

    return output;
}

/* Shifts the string along given a new 3' base. Pops off the base at the 5'
 * end, allowing for efficient scanning of a string for our sequence.
 */

constexpr seqhash::word LEAST_SIG=(1 << BITS_PER_BASE) - 1;

void shift_sequence(seqhash& in, const size_t n, char new3) {
    if (in.size()==0) {
        return;
    }

    size_t x=0;
    auto current=in[x];
    current >>= BITS_PER_BASE; // Remove least significant bits.

    while (x+1 < in.size()) {
        auto next=in[x+1];
        auto popped=next & LEAST_SIG; // Get least significant 2 bits.
        next >>= BITS_PER_BASE; // Remove least significant 2 bits.

        // Add 2 most significant bits.
        popped <<= (WIDTH_IN_BASES-1)*BITS_PER_BASE;
        current |= popped;
        in[x] = current;

        current=next;
        ++x;
    }

    // Adding 2 most significant bits (based on 'n', not the integer width).
    current |= (bitify(new3) << ((n - 1) % WIDTH_IN_BASES) * BITS_PER_BASE);
    in[x] = current;

    return;
}
