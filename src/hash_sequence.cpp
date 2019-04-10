#include "hash_sequence.h"

#ifdef DEBUG
#include <iostream>

void disgorge(const std::u32string& in) {
    for (size_t i=0; i<in.size(); ++i) {
        std::cout << static_cast<uint32_t>(in[i]) << ", ";
    }
    std::cout << std::endl;
    return;
}

#endif

int bitify (char base) {
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

constexpr int WIDTH_IN_BASES=16; // 16 bases, 2 bits each => char32_t.
constexpr int BITS_PER_BASE=2;

std::u32string hash_sequence(const char* p, const size_t n) { 
    int effective_length=(n + WIDTH_IN_BASES - 1)/WIDTH_IN_BASES;
    std::u32string output;
    output.reserve(effective_length);

    size_t i=0;
    while (i < n) {
        size_t limit=std::min(n, i+WIDTH_IN_BASES);

        uint32_t current=0;
        int shift=0;
        for (; i<limit; ++i) {
            current+=(bitify(p[i]) << shift);
            shift+=BITS_PER_BASE;
        }

        output.push_back(static_cast<char32_t>(current));
    }

    return output;
}

/* Shifts the string along given a new 3' base.
 * Pops off the base at the 5' end, allowing for
 * efficient scanning of a string for our sequence.
 */

void shift_sequence(std::u32string& in, const size_t n, char new3) {
    if (in.size()==0) {
        throw std::runtime_error("empty string cannot be shifted");
    }

    size_t x=0;
    uint32_t current=static_cast<uint32_t>(in[x]);
    current >>= BITS_PER_BASE; // Remove least significant bits.

    while (x+1 < in.size()) {
        uint32_t next=static_cast<uint32_t>(in[x+1]);
        uint32_t popped=next & 0x3; // Get least significant 2 bits.
        next >>= BITS_PER_BASE; // Remove least significant 2 bits.

        // Add 2 most significant bits.
        popped <<= (WIDTH_IN_BASES-1)*BITS_PER_BASE;
        current += popped; 
        in[x] = static_cast<char32_t>(current);

        current=next;
        ++x;
    }

    // Adding 2 most significant bits (based on 'n', not the integer width).
    current += (bitify(new3) << ((n - 1) % WIDTH_IN_BASES) * BITS_PER_BASE);
    in[x] = static_cast<char32_t>(current);

    return;
}

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
    ++ptr;
    const char next=*(ptr+len);
    shift_sequence(hashed, len, next);
    nvalid+=is_valid(next);
    return;
}

const std::u32string& hash_scanner::hash() const { return hashed; }

bool hash_scanner::valid() const { return nvalid==len; }

