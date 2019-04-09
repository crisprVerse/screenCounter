#include "hash_string.h"

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

/* Shifts the string along given a new 5' base.
 * Pops off the base at the 3' end.
 */

void shift_string(std::u32string& in, const size_t n, char new5) {
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
    current += (bitify(new5) << ((n - 1) % WIDTH_IN_BASES) * BITS_PER_BASE);
    in[x] = static_cast<char32_t>(current);

    return;
}

