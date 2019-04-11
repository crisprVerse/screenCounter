#include "hash_sequence.h"
#include <stdexcept>

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

constexpr uint32_t LEAST_SIG=(1 << BITS_PER_BASE) - 1;

void shift_sequence(std::u32string& in, const size_t n, char new3) {
    if (in.size()==0) {
        throw std::runtime_error("empty string cannot be shifted");
    }

    size_t x=0;
    uint32_t current=static_cast<uint32_t>(in[x]);
    current >>= BITS_PER_BASE; // Remove least significant bits.

    while (x+1 < in.size()) {
        uint32_t next=static_cast<uint32_t>(in[x+1]);
        uint32_t popped=next & LEAST_SIG; // Get least significant 2 bits.
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

/* Substitutes a base in the current hash string. 
 * This is done in a rolling manner, where the next base
 * is substituted upon the user calling advance().
 */ 

rolling_substitution::rolling_substitution(const std::u32string& in, size_t n) :
    hash(in), len(n), word(0), pos(0), state(0), current(0), original(0)
{
    if (len==0) { 
        throw std::runtime_error("hashed string must have positive length");
    }

    // Starting off by replacing the first base with 'A'.
    current=static_cast<uint32_t>(in[word]);
    original=current & LEAST_SIG;
    if (original==0) {
        advance();
    } else {
        current &= ~LEAST_SIG;
        hash[word]=static_cast<char32_t>(current);
    }
    return; 
}

bool rolling_substitution::advance() {
    ++state;
    if (state==4) {
        // Restoring original if we're moving along.
        current &= ~(LEAST_SIG << pos * BITS_PER_BASE);
        current |= (original << pos * BITS_PER_BASE);
        state=0;
        ++pos;

        if (pos+WIDTH_IN_BASES*word==len) {
            return false;
        } else if (pos==WIDTH_IN_BASES) {
            pos=0;
            hash[word]=static_cast<char32_t>(current);
            ++word;
            current=static_cast<uint32_t>(hash[word]);
        }

        // Updating original.
        original=current & (LEAST_SIG << pos * BITS_PER_BASE);
        original >>= pos * BITS_PER_BASE;
    } 

    if (state==original) {
        return advance();
    } 
    
    // Replacing with next state.
    current &= ~(LEAST_SIG << pos * BITS_PER_BASE);
    current |= (state << pos * BITS_PER_BASE);
    hash[word]=static_cast<char32_t>(current);
    return true;
}

const std::u32string& rolling_substitution::get() const { return hash; }

/* Deletes a base in the current hash string. 
 * This starts from the last base and rolls towards
 * the first base (simply because it's more convenient).
 */ 

rolling_deletion::rolling_deletion(const std::u32string& in, size_t n) :
    hash(in), word((n-1)/WIDTH_IN_BASES), pos((n-1)%WIDTH_IN_BASES), current(0), discarded(0)
{
    if (n==0) { 
        throw std::runtime_error("hashed string must have positive length");
    }

    // Starting off by deleting the last base.
    current=static_cast<uint32_t>(in[word]);
    discarded=current & (LEAST_SIG << pos * BITS_PER_BASE);
    current &= ~(LEAST_SIG << pos * BITS_PER_BASE);

    if (pos==0) {
        hash.pop_back();
    } else {
        hash[word]=static_cast<char32_t>(current);
    }
    return; 
}

bool rolling_deletion::advance() {
    if (pos==0) {
        if (word==0) {
            return false;
        }

        discarded <<= (WIDTH_IN_BASES - 1) * BITS_PER_BASE;
        pos=WIDTH_IN_BASES-1;
        --word;
        current=static_cast<uint32_t>(hash[word]);
    } else {
        --pos;
        discarded >>= BITS_PER_BASE;
    }

    // Restoring the discarded at its original position minus one.
    uint32_t tmp=current & (LEAST_SIG << pos * BITS_PER_BASE);
    current &= ~(LEAST_SIG << pos * BITS_PER_BASE);
    current |= discarded;
    discarded=tmp;
    hash[word]=static_cast<char32_t>(current);
    return true;
}

const std::u32string& rolling_deletion::get() const { return hash; }
