#ifndef KAORI_UTILS_HPP
#define KAORI_UTILS_HPP

#include <bitset>
#include <vector>
#include <array>

/**
 * @file utils.hpp
 * @brief Utilites for sequence matching.
 */

namespace kaori {

/**
 * How to deal with duplicated sequences in the pool of known barcodes.
 *
 * - `FIRST`: keep the first encountered instance of a set of duplicates.
 * - `LAST`: keep the last encountered instance of a set of duplicates.
 * - `NONE`: do not keep any duplicate sequences.
 * - `ERROR`: throw an error upon observing a duplicate sequence.
 */
enum class DuplicateAction : char { FIRST, LAST, NONE, ERROR };

/**
 * Strands of the read sequence to search.
 */
enum class SearchStrand : char { FORWARD, REVERSE, BOTH };

/**
 * @cond
 */
inline bool search_forward(SearchStrand x) {
    return x == SearchStrand::FORWARD || x == SearchStrand::BOTH;
}

inline bool search_reverse(SearchStrand x) {
    return x == SearchStrand::REVERSE || x == SearchStrand::BOTH;
}

template<bool allow_n_ = false, bool allow_iupac_ = false>
char complement_base(char b) {
    char output;
    switch (b) {
        case 'A': case 'a':
            output = 'T';
            break;
        case 'C': case 'c':
            output = 'G';
            break;
        case 'G': case 'g':
            output = 'C';
            break;
        case 'T': case 't':
            output = 'A';
            break;

        case 'N': case 'n':
            if constexpr(allow_n_ || allow_iupac_) {
                output = 'N';
                break;
            }

        case 'R': case 'r':
            if constexpr(allow_iupac_) {
                output = 'Y';
                break;
            }
        case 'Y': case 'y':
            if constexpr(allow_iupac_) {
                output = 'R';
                break;
            }
        case 'S': case 's':
            if constexpr(allow_iupac_) {
                output = 'S'; // S = A/T, so complement is just itself.
                break;
            }
        case 'W': case 'w':
            if constexpr(allow_iupac_) {
                output = 'W'; // W = C/G, so complement is just itself.
                break;
            }
        case 'K': case 'k':
            if constexpr(allow_iupac_) {
                output = 'M';
                break;
            }
        case 'M': case 'm':
            if constexpr(allow_iupac_) {
                output = 'K';
                break;
            }

        case 'B': case 'b':
            if constexpr(allow_iupac_) {
                output = 'V'; // B can't be A, so complement can't be T ==> V.
                break;
            }
        case 'D': case 'd':
            if constexpr(allow_iupac_) {
                output = 'H'; // D can't be C, so complement can't be G ==> H.
                break;
            }
        case 'H': case 'h':
            if constexpr(allow_iupac_) {
                output = 'D'; // H can't be G, so complement can't be C ==> D.
                break;
            }
        case 'V': case 'v':
            if constexpr(allow_iupac_) {
                output = 'B'; // V can't be T, so complement can't be A ==> B.
                break;
            }

        default:
            throw std::runtime_error("cannot complement unknown base '" + std::string(1, b) + "'");
    }
    return output;
}

inline bool is_standard_base(char b) {
    bool okay = false;
    switch (b) {
        case 'A': case 'a': 
        case 'C': case 'c': 
        case 'G': case 'g': 
        case 'T': case 't':
            okay = true;
            break;
    }
    return okay;
}

template<size_t N>
void shift_hash(std::bitset<N>& x) {
    x <<= 4;
}

template<size_t N>
void add_base_to_hash(std::bitset<N>& x, char b) {
    shift_hash(x);
    switch (b) {
        case 'A': case 'a':
            x.set(0);
            break;
        case 'C': case 'c':
            x.set(1);
            break;
        case 'G': case 'g':
            x.set(2);
            break;
        case 'T': case 't':
            x.set(3);
            break;
        default:
            throw std::runtime_error("unknown base '" + std::string(1, b) + "'");
            break;
    }
    return;
}

template<size_t N>
void add_other_to_hash(std::bitset<N>& x) {
    shift_hash(x);
    x.set(0);
    x.set(1);
    x.set(2);
    x.set(3);
    return;
}

template<size_t V>
void sort_combinations(std::vector<std::array<int, V> >& combinations, const std::array<size_t, V>& num_options) {
    // Going back to front as the last iteration gives the slowest changing index.
    // This ensures that we get the same results as std::sort() on the arrays.
    for (size_t i_ = 0; i_ < V; ++i_) {
        auto i = V - i_ - 1;

        std::vector<size_t> counts(num_options[i] + 1);
        for (const auto& x : combinations) {
            ++(counts[x[i] + 1]);
        }

        for (size_t j = 1; j < counts.size(); ++j) {
            counts[j] += counts[j-1];
        }

        std::vector<std::array<int, V> > copy(combinations.size());
        for (const auto& x : combinations) {
            auto& pos = counts[x[i]];
            copy[pos] = x;
            ++pos;
        }

        combinations.swap(copy);
    }
}
/**
 * @endcond
 */

}

#endif
