#ifndef KAORI_MISMATCH_TRIE_HPP
#define KAORI_MISMATCH_TRIE_HPP

#include <vector>
#include <stdexcept>
#include "utils.hpp"

namespace kaori {

class MismatchTrie {
public:
    MismatchTrie(size_t n = 0) : length(n), pointers(4, -1), counter(0) {}

    MismatchTrie(const std::vector<const char*>& seq, size_t n) : MismatchTrie(n) {
        for (auto s : seq) {
            add(s);
        }
    }

public:
    void add(const char* seq) {
        int position = 0;

        for (size_t i = 0; i < length; ++i) {
            auto& current = pointers[position + base_shift(seq[i])];

            if (i + 1 == length) {
                // Last position is the index of the sequence.
                if (current >= 0) {
                    throw std::runtime_error("duplicate sequences detected when constructing the trie");
                }
                current = counter;
            } else {
                if (current < 0) {
                    current = pointers.size();
                    position = current;
                    pointers.resize(position + 4, -1);
                } else {
                    position = current;
                }
            }
        }

        ++counter;
    }

public:
    std::pair<int, int> search(const char* seq, int max_mismatch) const {
        return search(seq, 0, 0, 0, max_mismatch);
    }

    std::pair<int, int> search(const char* seq, size_t pos, int node, int mismatches, int& max_mismatch) const {
        int shift = base_shift(seq[pos]);
        int current = pointers[node + shift];

        // At the end: we prepare to return the actual values. We also refine
        // the max number of mismatches so that we don't search for things with
        // more mismatches than the best hit that was already encountered.
        if (pos + 1 == length) {
            if (current >= 0) {
                max_mismatch = mismatches;
                return std::make_pair(current, mismatches);
            }

            int alt = -1;
            ++mismatches;
            if (mismatches <= max_mismatch) {
                bool found = false;
                for (int s = 0; s < 4; ++s) {
                    if (shift == s) { 
                        continue;
                    }

                    int candidate = pointers[node + s];
                    if (candidate >= 0) {
                        if (found) { // ambiguous, so we quit early.
                            alt = -1;
                            break;
                        }
                        alt = candidate;
                        max_mismatch = mismatches;
                        found = true;
                    }
                }
            }
            return std::make_pair(alt, mismatches);

        } else {
            ++pos;

            std::pair<int, int> best(-1, max_mismatch + 1);
            if (current >= 0) {
                best = search(seq, pos, current, mismatches, max_mismatch);
            }

            ++mismatches;
            if (mismatches <= max_mismatch) {
                bool found = false;
                for (int s = 0; s < 4; ++s) {
                    if (shift == s) { 
                        continue;
                    } 
                    
                    int alt = pointers[node + s];
                    if (alt < 0) {
                        continue;
                    }

                    auto chosen = search(seq, pos, alt, mismatches, max_mismatch);
                    if (chosen.second < best.second) {
                        best = chosen;
                    } else if (chosen.second == best.second) {
                        best.first = -1;
                    }
                }
            }

            return best;
        }
    }

private:
    size_t length;
    std::vector<int> pointers;
    int counter;

    static int base_shift(char base) {
        int shift = 0;
        switch (base) {
            case 'A': case 'a':
                break;
            case 'C': case 'c':
                shift = 1;
                break;
            case 'G': case 'g':
                shift = 2;
                break;
            case 'T': case 't':
                shift = 3;
                break;
            default:
                throw std::runtime_error("unknown base '" + std::string(1, base) + "' detected when constructing the trie");
        }
        return shift;
    }
};

}

#endif
