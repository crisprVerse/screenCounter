#ifndef KAORI_MISMATCH_TRIE_HPP
#define KAORI_MISMATCH_TRIE_HPP

#include <array>
#include <vector>
#include <stdexcept>
#include <numeric>
#include "utils.hpp"

namespace kaori {

class MismatchTrie {
public:
    MismatchTrie(size_t n = 0) : length(n), pointers(4, -1), counter(0) {}

    MismatchTrie(const std::vector<const char*>& seq, size_t n, bool duplicates = false) : MismatchTrie(n) {
        for (auto s : seq) {
            add(s, duplicates);
        }
    }

public:
    void add(const char* seq, bool duplicates = false) {
        int position = 0;

        for (size_t i = 0; i < length; ++i) {
            auto& current = pointers[position + base_shift(seq[i])];

            if (i + 1 == length) {
                // Last position is the index of the sequence.
                if (current >= 0) {
                    if (!duplicates) {
                        throw std::runtime_error("duplicate sequences detected when constructing the trie");
                    }
                } else {
                    current = counter;
                }
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

    size_t get_length() const {
        return length;
    }

protected:
    size_t length;
    std::vector<int> pointers;

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

private:
    int counter;
};

class SimpleMismatchTrie : public MismatchTrie {
public:
    SimpleMismatchTrie(size_t n = 0) : MismatchTrie(n) {}

    SimpleMismatchTrie(const std::vector<const char*>& seq, size_t n, bool duplicates = false) : MismatchTrie(seq, n, duplicates) {}

    std::pair<int, int> search(const char* seq, int max_mismatch) const {
        return search(seq, 0, 0, 0, max_mismatch);
    }
private:
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
};

template<size_t num_segments>
class SegmentedMismatchTrie : public MismatchTrie {
public:
    SegmentedMismatchTrie() {}

    SegmentedMismatchTrie(std::array<int, num_segments> segments) : MismatchTrie(std::accumulate(segments.begin(), segments.end(), 0)), boundaries(segments) {
        for (size_t i = 1; i < num_segments; ++i) {
            boundaries[i] += boundaries[i-1];
        }
    }

    SegmentedMismatchTrie(const std::vector<const char*>& seq, std::array<int, num_segments> segments, bool duplicates = false) : SegmentedMismatchTrie(segments) {
        for (auto s : seq) {
            add(s, duplicates);
        }
    }

    struct SearchResult {
        SearchResult() : per_segment() {}
        int index = 0;
        int total = 0;
        std::array<int, num_segments> per_segment;
    };

    SearchResult search(const char* seq, const std::array<int, num_segments>& mismatches) const {
        int max_mismatches = std::accumulate(mismatches.begin(), mismatches.end(), 0);
        return search(seq, 0, 0, SearchResult(), mismatches, max_mismatches);
    }

private:
    SearchResult search(
        const char* seq, 
        size_t pos, 
        size_t segment_id,
        SearchResult state,
        const std::array<int, num_segments>& segment_mismatches, 
        int& max_mismatch
    ) const {
        // Note that, during recursion, state.index does double duty 
        // as the index of the node on the trie.
        int node = state.index;

        int shift = base_shift(seq[pos]);
        int current = pointers[node + shift];

        // At the end: we prepare to return the actual values. We also refine
        // the max number of mismatches so that we don't search for things with
        // more mismatches than the best hit that was already encountered.
        if (pos + 1 == length) {
            if (current >= 0) {
                max_mismatch = state.total;
                state.index = current;
                return state;
            }

            state.index = -1;
            ++state.total;
            auto& current_segment_mm = state.per_segment[segment_id];
            ++current_segment_mm;

            if (state.total <= max_mismatch && current_segment_mm <= segment_mismatches[segment_id]) {
                bool found = false;
                for (int s = 0; s < 4; ++s) {
                    if (shift == s) { 
                        continue;
                    }

                    int candidate = pointers[node + s];
                    if (candidate >= 0) {
                        if (found) { // ambiguous, so we quit early.
                            state.index = -1;
                            break;
                        }
                        state.index = candidate;
                        max_mismatch = state.total;
                        found = true;
                    }
                }
            }
            return state;

        } else {
            ++pos;
            if (pos == boundaries[segment_id]) {
                ++segment_id;
            }

            SearchResult best;
            best.index = -1;
            best.total = max_mismatch + 1;

            if (current >= 0) {
                state.index = current;
                best = search(seq, pos, segment_id, state, segment_mismatches, max_mismatch);
            }

            ++state.total;
            auto& current_segment_mm = state.per_segment[segment_id];
            ++current_segment_mm;

            if (state.total <= max_mismatch && current_segment_mm <= segment_mismatches[segment_id]) {
                bool found = false;
                for (int s = 0; s < 4; ++s) {
                    if (shift == s) { 
                        continue;
                    } 
                    
                    int alt = pointers[node + s];
                    if (alt < 0) {
                        continue;
                    }

                    state.index = alt;
                    auto chosen = search(seq, pos, segment_id, state, segment_mismatches, max_mismatch);
                    if (chosen.total < best.total) {
                        best = chosen;
                    } else if (chosen.total == best.total) { // ambiguous
                        best.index = -1;
                    }
                }
            }

            return best;
        }
    }
private:
    std::array<int, num_segments> boundaries;
};

}

#endif
