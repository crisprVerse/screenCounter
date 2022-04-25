#ifndef KAORI_MISMATCH_TRIE_HPP
#define KAORI_MISMATCH_TRIE_HPP

#include <array>
#include <vector>
#include <stdexcept>
#include <numeric>
#include "utils.hpp"
#include "BarcodePool.hpp"

/**
 * @file MismatchTrie.hpp
 *
 * @brief Defines the `MismatchTrie` class and its subclasses.
 */

namespace kaori {

/**
 * @brief Base class for the mismatch search.
 *
 * Given a (typically read-derived) sequence, we can perform a mismatch-aware search to known sequences in a pool of barcode sequences.
 * The idea is to find the barcode with the fewest mismatches to the input sequence.
 * Any number of mismatches are supported; subclasses will decide how the mismatches can be distributed throughout the length of the sequence.
 */
class MismatchTrie {
public:
    /**
     * @param barcode_length Length of the barcodes in the pool.
     */
    MismatchTrie(size_t barcode_length = 0) : length(barcode_length), pointers(4, -1), counter(0) {}

    /**
     * @param barcode_pool Pool of known barcode sequences.
     * @param duplicates Whether duplicated sequences in `barcode_pool` should be supported, see `add()`.
     */
    MismatchTrie(const BarcodePool& barcode_pool, bool duplicates = false) : MismatchTrie(barcode_pool.length) {
        for (auto s : barcode_pool.pool) {
            add(s, duplicates);
        }
    }

public:
    /**
     * @param[in] barcode_seq Pointer to a character array containing a barcode sequence.
     * The array should have length equal to `get_length()`.
     * @param duplicates Whether duplicate sequences are allowed.
     * If `false`, an error is raised if `seq` is a duplicate of a previously `add()`ed sequence.
     * If `true`, only the first instance of the duplicates will be reported in searches.
     *
     * @return The barcode sequence is added to the trie.
     * The index of the newly added sequence is defined as the number of sequences that were previously added. 
     */
    void add(const char* barcode_seq, bool duplicates = false) {
        int position = 0;

        for (size_t i = 0; i < length; ++i) {
            auto& current = pointers[position + base_shift(barcode_seq[i])];

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

    /**
     * @return The length of the barcode sequences.
     */
    size_t get_length() const {
        return length;
    }

    /**
     * @return The number of barcode sequences added.
     */
    int size() const {
        return counter;
    }

protected:
    /**
     * @cond
     */
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
    /**
     * @endcond
     */

private:
    int counter;
};

/**
 * @brief Search for barcodes with mismatches anywhere.
 *
 * This `MismatchTrie` subclass will search for the best match to known sequences in a barcode pool.
 * Any number of mismatches are supported, distributed anywhere throughout the sequence. 
 */
class AnyMismatches : public MismatchTrie {
public:
    /**
     * @param barcode_length Length of the barcode sequences.
     */
    AnyMismatches(size_t barcode_length = 0) : MismatchTrie(barcode_length) {}

    /**
     * @param barcode_pool Pool of known barcode sequences.
     * @param duplicates Whether duplicated sequences in `barcode_pool` should be supported, see `add()`.
     */
    AnyMismatches(const BarcodePool& barcode_pool, bool duplicates = false) : MismatchTrie(barcode_pool, duplicates) {}

public:
    /**
     * @param[in] search_seq Pointer to a character array containing a sequence to use for searching the barcode pool.
     * This is assumed to be of length equal to `get_length()` and is typically derived from a read.
     * @param max_mismatches Maximum number of mismatches in the search.
     *
     * @return Pair containing:
     * 1. The index of the barcode sequence with the lowest number of mismatches to `search_seq`.
     *    If multiple sequences have the same lowest number of mismatches, the match is ambiguous and -1 is returned.
     *    If all sequences have more mismatches than `max_mismatches`, -1 is returned.
     * 2. The number of mismatches.
     */
    std::pair<int, int> search(const char* search_seq, int max_mismatches) const {
        return search(search_seq, 0, 0, 0, max_mismatches);
    }

private:
    std::pair<int, int> search(const char* seq, size_t pos, int node, int mismatches, int& max_mismatches) const {
        int shift = base_shift(seq[pos]);
        int current = pointers[node + shift];

        // At the end: we prepare to return the actual values. We also refine
        // the max number of mismatches so that we don't search for things with
        // more mismatches than the best hit that was already encountered.
        if (pos + 1 == length) {
            if (current >= 0) {
                max_mismatches = mismatches;
                return std::make_pair(current, mismatches);
            }

            int alt = -1;
            ++mismatches;
            if (mismatches <= max_mismatches) {
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
                        max_mismatches = mismatches;
                        found = true;
                    }
                }
            }
            return std::make_pair(alt, mismatches);

        } else {
            ++pos;

            std::pair<int, int> best(-1, max_mismatches + 1);
            if (current >= 0) {
                best = search(seq, pos, current, mismatches, max_mismatches);
            }

            ++mismatches;
            if (mismatches <= max_mismatches) {
                bool found = false;
                for (int s = 0; s < 4; ++s) {
                    if (shift == s) { 
                        continue;
                    } 
                    
                    int alt = pointers[node + s];
                    if (alt < 0) {
                        continue;
                    }

                    auto chosen = search(seq, pos, alt, mismatches, max_mismatches);
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

/**
 * @brief Search for barcodes with segmented mismatches.
 *
 * This `MismatchTrie` subclass will search for the best match to known sequences in a barcode pool.
 * However, the distribution of mismatches is restricted in different segments of the sequence, e.g., 1 mismatch in the first 4 bp, 3 mismatches for the next 10 bp, and so on.
 * The intention is to enable searching for concatenations of variable region sequences (and barcodes), where each segment is subject to a different number of mismatches.
 *
 * @tparam num_segments Number of segments to consider.
 */
template<size_t num_segments>
class SegmentedMismatches : public MismatchTrie {
public:
    /**
     * Default constructor.
     */
    SegmentedMismatches() {}

    /**
     * @param segments Length of each segment of the sequence.
     * Each entry should be positive and the sum should be equal to the total length of the barcode sequence.
     */
    SegmentedMismatches(std::array<int, num_segments> segments) : MismatchTrie(std::accumulate(segments.begin(), segments.end(), 0)), boundaries(segments) {
        for (size_t i = 1; i < num_segments; ++i) {
            boundaries[i] += boundaries[i-1];
        }
    }

    /**
     * @param barcode_pool Possible set of known sequences for the variable region.
     * @param segments Length of each segment of the sequence.
     * Each entry should be positive and the sum should be equal to the total length of the barcode sequence.
     * @param duplicates Whether duplicated sequences in `barcode_pool` should be supported, see `add()`.
     */
    SegmentedMismatches(const BarcodePool& barcode_pool, std::array<int, num_segments> segments, bool duplicates = false) : SegmentedMismatches(segments) {
        if (length != barcode_pool.length) {
            throw std::runtime_error("length of barcode sequences should equal the sum of segment lengths");
        }
        for (auto s : barcode_pool.pool) {
            add(s, duplicates);
        }
    }

public:
    /**
     * @brief Result of the segmented search.
     */
    struct Result {
        /**
         * @cond
         */
        Result() : per_segment() {}
        /**
         * @endcond
         */

        /**
         * Index of the known barcode sequence matching the input sequence in `search()`.
         * If no unambiguous match is found, -1 is reported.
         */
        int index = 0;

        /**
         * Total number of mismatches between the barcode sequence from `index` and the input sequence.
         */
        int total = 0;

        /**
         * Number of mismatches in each segment of the sequence.
         */
        std::array<int, num_segments> per_segment;
    };

    /**
     * @param[in] search_seq Pointer to a character array containing a sequence to use for searching the barcode pool.
     * This is assumed to be of length equal to `get_length()` and is typically derived from a read.
     * @param max_mismatches Maximum number of mismatches for each segment.
     * Each entry should be non-negative.
     *
     * @return A `Result` containing the index of the barcode sequence where the number of mismatches in each segment is less than or equal to `max_mismatches`.
     * - If multiple barcode sequences satisfy this condition, the barcode sequence with the lowest total number of mismatches is reported.
     * - If multiple barcode sequences share the same lowest total, the match is ambiguous and -1 is reported.
     * - If no barcode sequences satisfy the `max_mismatches` condition, -1 is reported.
     */
    Result search(const char* search_seq, const std::array<int, num_segments>& max_mismatches) const {
        int total_mismatches = std::accumulate(max_mismatches.begin(), max_mismatches.end(), 0);
        return search(search_seq, 0, 0, Result(), max_mismatches, total_mismatches);
    }

private:
    Result search(
        const char* seq, 
        size_t pos, 
        size_t segment_id,
        Result state,
        const std::array<int, num_segments>& segment_mismatches, 
        int& total_mismatches
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
                total_mismatches = state.total;
                state.index = current;
                return state;
            }

            state.index = -1;
            ++state.total;
            auto& current_segment_mm = state.per_segment[segment_id];
            ++current_segment_mm;

            if (state.total <= total_mismatches && current_segment_mm <= segment_mismatches[segment_id]) {
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
                        total_mismatches = state.total;
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

            Result best;
            best.index = -1;
            best.total = total_mismatches + 1;

            if (current >= 0) {
                state.index = current;
                best = search(seq, pos, segment_id, state, segment_mismatches, total_mismatches);
            }

            ++state.total;
            auto& current_segment_mm = state.per_segment[segment_id];
            ++current_segment_mm;

            if (state.total <= total_mismatches && current_segment_mm <= segment_mismatches[segment_id]) {
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
                    auto chosen = search(seq, pos, segment_id, state, segment_mismatches, total_mismatches);
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
