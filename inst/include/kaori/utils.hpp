#ifndef KAORI_UTILS_HPP
#define KAORI_UTILS_HPP

#include <bitset>
#include <vector>
#include <array>

namespace {

template<size_t N>
constexpr std::bitset<N> A_(1);

template<size_t N>
constexpr std::bitset<N> C_(2);

template<size_t N>
constexpr std::bitset<N> G_(4);

template<size_t N>
constexpr std::bitset<N> T_(8);

template<size_t N>
constexpr std::bitset<N> other_(15);

inline char reverse_complement(char b) {
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
        default:
            throw std::runtime_error("unknown base '" + std::string(1, b) + "'");
            break;
    }
    return output;
}

template<size_t N>
void shift(std::bitset<N>& x) {
    x <<= 4;
}

inline bool is_good(char b) {
    bool okay = false;
    switch (b) {
        case 'A': case 'a':
            okay = true;
            break;
        case 'C': case 'c':
            okay = true;
            break;
        case 'G': case 'g':
            okay = true;
            break;
        case 'T': case 't':
            okay = true;
            break;
    }
    return okay;
}

template<size_t N>
void add_base(std::bitset<N>& x, char b) {
    shift(x);
    switch (b) {
        case 'A': case 'a':
            x |= A_<N>;
            break;
        case 'C': case 'c':
            x |= C_<N>;
            break;
        case 'G': case 'g':
            x |= G_<N>;
            break;
        case 'T': case 't':
            x |= T_<N>;
            break;
        default:
            throw std::runtime_error("unknown base '" + std::string(1, b) + "'");
            break;
    }
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

}

#endif
