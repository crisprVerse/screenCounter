#ifndef KAORI_SIMPLE_SINGLE_MATCH_HPP
#define KAORI_SIMPLE_SINGLE_MATCH_HPP

#include "ScanTemplate.hpp"
#include "BarcodePool.hpp"
#include "BarcodeSearch.hpp"
#include "utils.hpp"

#include <string>
#include <unordered_map>
#include <vector>

/**
 * @file SimpleSingleMatch.hpp
 *
 * @brief Defines the `SimpleSingleMatch` class.
 */

namespace kaori {

/**
 * @brief Search for a template with a single variable region.
 *
 * This class implements the most common use case for barcode matching, where the template sequence has a single variable region.
 * It will find a match to any valid barcoding element, i.e., the realization of the template where the variable region is replaced with one sequence from a pool of known barcodes.
 * No restrictions are placed on the distribution of mismatches throughout the barcoding element.
 * 
 * @tparam max_size Maximum length of the template sequence.
 */
template<size_t max_size>
class SimpleSingleMatch {
public:
    /**
     * @brief Optional parameters for `SimpleSingleMatch`.
     */
    struct Options {
        /**
         * Maximum number of mismatches for any search performed by `SimpleBarcodeSearch::search`.
         */
        int max_mismatches = 0;

        /** 
         * How duplicated barcode sequences should be handled.
         */
        DuplicateAction duplicates = DuplicateAction::ERROR;

        /** 
         * Strand(s) of the read sequence to search.
         */
        SearchStrand strand = SearchStrand::FORWARD;
    };

public:
    /**
     * @param[in] template_seq Pointer to a character array containing the template sequence, see `ScanTemplate`.
     * @param template_length Length of the array pointed to by `barcode_length`.
     * This should be less than or equal to `max_size`.
     * @param barcode_pool Known sequences for the single variable region in `template_seq`.
     * @param options Optional parameters.
     */
    SimpleSingleMatch(
        const char* template_seq, 
        size_t template_length, 
        const BarcodePool& barcode_pool, 
        const Options& options
    ) : 
        num_options(barcode_pool.pool.size()),
        forward(search_forward(options.strand)), 
        reverse(search_reverse(options.strand)),
        max_mm(options.max_mismatches),
        constant(template_seq, template_length, options.strand)
    {
        // Exact strandedness doesn't matter here, just need the number and length.
        const auto& regions = constant.variable_regions();
        if (regions.size() != 1) {
            throw std::runtime_error("expected one variable region in the constant template");
        }

        size_t var_length = regions[0].second - regions[0].first;
        if (var_length != barcode_pool.length) {
            throw std::runtime_error("length of barcode_pool sequences (" + std::to_string(barcode_pool.length) + 
                ") should be the same as the barcode_pool region (" + std::to_string(var_length) + ")");
        }

        SimpleBarcodeSearch::Options bs_opt;
        bs_opt.duplicates = options.duplicates;
        bs_opt.max_mismatches = max_mm;

        if (forward) {
            bs_opt.reverse = false;
            forward_lib = SimpleBarcodeSearch(barcode_pool, bs_opt);
        }
        if (reverse) {
            bs_opt.reverse = true;
            reverse_lib = SimpleBarcodeSearch(barcode_pool, bs_opt);
        }
    }

public:
    /**
     * @brief State of the search on a read sequence.
     */
    struct State {
        /**
         * Index of the known barcode sequence that matches the variable region in the read sequence.
         * This will be -1 if no match was found.
         */
        int index = 0;

        /**
         * Position of the match to the template, reported as the position on the read at the start of the template.
         * This should only be used if `index != -1`.
         */
        size_t position = 0;

        /**
         * Total number of mismatches, including both the constant and variable regions.
         * This should only be used if `index != -1`.
         */
        int mismatches = 0;

        /**
         * Total number of mismatches in the variable region.
         * This should only be used if `index != -1`.
         */
        int variable_mismatches = 0;

        /**
         * Whether the match was found on the reverse strand of the read sequence.
         */
        bool reverse = false;

        /**
         * @cond
         */
        typename SimpleBarcodeSearch::State forward_details, reverse_details;
        /**
         * @endcond
         */
    };

    /**
     * Initialize the search state for thread-safe execution.
     *
     * @return A new `State` object.
     */
    State initialize() const {
        return State();
    }

    /**
     * Incorporate search cache optimizations from `state`, see `SimpleBarcodeSearch::reduce()` for details.
     * This allows regular consolidation and sharing of optimizations across threads.
     * On return, optimizations from `state` are incorporated into this `SimpleSingleMatch` instance.
     *
     * @param state A state object generated by `initialize()`.
     * Typically this has been used in `search_first()` or `search_best()` at least once.
     */
    void reduce(State& state) {
        if (forward) {
            forward_lib.reduce(state.forward_details);
        }
        if (reverse) {
            reverse_lib.reduce(state.reverse_details);
        }
    }

private:
    bool has_match(int obs_mismatches) const {
        return (obs_mismatches >= 0 && obs_mismatches <= max_mm);
    }

    void forward_match(const char* seq, const typename ScanTemplate<max_size>::State& details, State& state) const {
        auto start = seq + details.position;
        const auto& range = constant.variable_regions()[0];
        std::string curseq(start + range.first, start + range.second);
        forward_lib.search(curseq, state.forward_details, max_mm - details.forward_mismatches);
    }

    void reverse_match(const char* seq, const typename ScanTemplate<max_size>::State& details, State& state) const {
        auto start = seq + details.position;
        const auto& range = constant.template variable_regions<true>()[0];
        std::string curseq(start + range.first, start + range.second);
        reverse_lib.search(curseq, state.reverse_details, max_mm - details.reverse_mismatches);
    }

public:
    /**
     * Search a read for the first match to a valid barcoding element.
     * A match is only reported if the number of mismatches of the entire barcoding element's sequence to the read is no greater than `max_mismatches` (see the constructor)
     * and there is exactly one barcode sequence with the fewest mismatches to the read sequence at the variable region.
     *
     * @param[in] read_seq Pointer to a character array containing the read sequence.
     * @param read_length Length of the read sequence.
     * @param state State object, used to store the search result.
     *
     * @return Whether an appropriate match was found.
     * If `true`, `state` is filled with the details of the first match.
     */
    bool search_first(const char* read_seq, size_t read_length, State& state) const {
        auto deets = constant.initialize(read_seq, read_length);
        bool found = false;
        state.index = -1;
        state.mismatches = 0;
        state.variable_mismatches = 0;

        auto update = [&](bool rev, int const_mismatches, const typename SimpleBarcodeSearch::State& x) -> bool {
            if (x.index < 0) {
                return false;
            }

            int total = const_mismatches + x.mismatches;
            if (total > max_mm) {
                return false;
            }

            found = true;
            state.position = deets.position;
            state.mismatches = total;
            state.reverse = rev;
            state.index = x.index;
            state.variable_mismatches = x.mismatches;
            return true;
        };

        while (!deets.finished) {
            constant.next(deets);

            if (forward && has_match(deets.forward_mismatches)) {
                forward_match(read_seq, deets, state);
                if (update(false, deets.forward_mismatches, state.forward_details)) {
                    break;
                }
            }

            if (reverse && has_match(deets.reverse_mismatches)) {
                reverse_match(read_seq, deets, state);
                if (update(true, deets.reverse_mismatches, state.reverse_details)) {
                    break;
                }
            }
        }

        return found;
    }

    /**
     * Search a read for the best match to a valid barcoding element. 
     * This is slower than `search_first()` but will find the matching position with the fewest mismatches.
     * If multiple positions are tied for the fewest mismatches, no match is reported.
     *
     * @param[in] read_seq Pointer to a character array containing the read sequence.
     * @param read_length Length of the read sequence.
     * @param state State object, used to store the search result.
     *
     * @return Whether a match was found.
     * If `true`, `state` is filled with the details of the best match.
     */
    bool search_best(const char* read_seq, size_t read_length, State& state) const {
        auto deets = constant.initialize(read_seq, read_length);
        state.index = -1;
        bool found = false;
        int best = max_mm + 1;

        auto update = [&](bool rev,  int const_mismatches, const typename SimpleBarcodeSearch::State& x) -> void {
            if (x.index < 0) {
                return;
            }

            auto total = x.mismatches + const_mismatches;
            if (total == best) { 
                if (state.index != x.index) { // ambiguous, setting back to a mismatch.
                    found = false;
                    state.index = -1;
                }
            } else if (total < best) {
                found = true;
                best = total; 
                // A further optimization at this point would be to narrow
                // max_mm to the current 'best'. But this probably
                // isn't worth it.

                state.index = x.index;
                state.mismatches = total;
                state.variable_mismatches = x.mismatches;
                state.position = deets.position;
                state.reverse = rev;
            }
        };

        while (!deets.finished) {
            constant.next(deets);

            if (forward && has_match(deets.forward_mismatches)) {
                forward_match(read_seq, deets, state);
                update(false, deets.forward_mismatches, state.forward_details);
            }

            if (reverse && has_match(deets.reverse_mismatches)) {
                reverse_match(read_seq, deets, state);
                update(true, deets.reverse_mismatches, state.reverse_details);
            }
        }

        return found;
    }

private:
    size_t num_options;
    bool forward, reverse;
    int max_mm;

    ScanTemplate<max_size> constant;
    SimpleBarcodeSearch forward_lib, reverse_lib;
};

}

#endif
