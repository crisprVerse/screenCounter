#ifndef KAORI_SIMPLE_SINGLE_MATCH_HPP
#define KAORI_SIMPLE_SINGLE_MATCH_HPP

#include "ConstantTemplate.hpp"
#include "VariableLibrary.hpp"
#include "utils.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace kaori {

template<size_t N>
class SimpleSingleMatch {
public:
    SimpleSingleMatch(const char* s, size_t n, bool f, bool r, const std::vector<const char*>& options, int mm = 0, bool duplicates = false) : 
        num_options(options.size()),
        forward(f), 
        reverse(r),
        max_mismatches(mm),
        constant(s, n, f, r)
    {
        // Exact strandedness doesn't matter here, just need the number and length.
        const auto& regions = constant.variable_regions();
        if (regions.size() != 1) {
            throw std::runtime_error("expected a single variable region only");
        }
        size_t var_length = regions[0].second - regions[0].first;

        if (forward) {
            forward_lib = SimpleVariableLibrary(options, var_length, max_mismatches, false, duplicates);
        }
        if (reverse) {
            reverse_lib = SimpleVariableLibrary(options, var_length, max_mismatches, true, duplicates);
        }
    }

public:
    struct SearchState {
        size_t position = 0;
        int index = 0;
        int mismatches = 0;
        int variable_mismatches = 0;
        bool reverse = false;

        /**
         * @cond
         */
        typename SimpleVariableLibrary::SearchState forward_details, reverse_details;
        /**
         * @endcond
         */
    };

    SearchState initialize() const {
        return SearchState();
    }

    void reduce(SearchState& state) {
        if (forward) {
            forward_lib.reduce(state.forward_details);
        }
        if (reverse) {
            reverse_lib.reduce(state.reverse_details);
        }
    }

private:
    bool has_match(int obs_mismatches) const {
        return (obs_mismatches >= 0 && obs_mismatches <= max_mismatches);
    }

    void forward_match(const char* seq, const typename ConstantTemplate<N>::MatchDetails& details, SearchState& state) const {
        auto start = seq + details.position;
        const auto& range = constant.variable_regions()[0];
        std::string curseq(start + range.first, start + range.second);
        forward_lib.match(curseq, state.forward_details, max_mismatches - details.forward_mismatches);
    }

    void reverse_match(const char* seq, const typename ConstantTemplate<N>::MatchDetails& details, SearchState& state) const {
        auto start = seq + details.position;
        const auto& range = constant.template variable_regions<true>()[0];
        std::string curseq(start + range.first, start + range.second);
        reverse_lib.match(curseq, state.reverse_details, max_mismatches - details.reverse_mismatches);
    }

public:
    bool search_first(const char* seq, size_t len, SearchState& state) const {
        auto deets = constant.initialize(seq, len);
        bool found = false;
        state.index = -1;
        state.mismatches = 0;
        state.variable_mismatches = 0;

        auto update = [&](bool rev, int const_mismatches, const typename SimpleVariableLibrary::SearchState& x) -> bool {
            if (x.index < 0) {
                return false;
            }

            int total = const_mismatches + x.mismatches;
            if (total > max_mismatches) {
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
                forward_match(seq, deets, state);
                if (update(false, deets.forward_mismatches, state.forward_details)) {
                    break;
                }
            }

            if (reverse && has_match(deets.reverse_mismatches)) {
                reverse_match(seq, deets, state);
                if (update(true, deets.reverse_mismatches, state.reverse_details)) {
                    break;
                }
            }
        }

        return found;
    }

    bool search_best(const char* seq, size_t len, SearchState& state) const {
        auto deets = constant.initialize(seq, len);
        state.index = -1;
        bool found = false;
        int best = max_mismatches + 1;

        auto update = [&](bool rev,  int const_mismatches, const typename SimpleVariableLibrary::SearchState& x) -> void {
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
                // max_mismatches to the current 'best'. But this probably
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
                forward_match(seq, deets, state);
                update(false, deets.forward_mismatches, state.forward_details);
            }

            if (reverse && has_match(deets.reverse_mismatches)) {
                reverse_match(seq, deets, state);
                update(true, deets.reverse_mismatches, state.reverse_details);
            }
        }

        return found;
    }

private:
    size_t num_options;
    bool forward, reverse;
    int max_mismatches;

    ConstantTemplate<N> constant;
    SimpleVariableLibrary forward_lib, reverse_lib;
};

}

#endif
