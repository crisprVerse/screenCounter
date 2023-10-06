#ifndef KAORI_COMBINATORIAL_BARCODES_SINGLE_END_HPP
#define KAORI_COMBINATORIAL_BARCODES_SINGLE_END_HPP

#include "../ScanTemplate.hpp"
#include "../BarcodeSearch.hpp"
#include "../utils.hpp"

#include <array>
#include <vector>

/**
 * @file CombinatorialBarcodesSingleEnd.hpp
 *
 * @brief Process single-end combinatorial barcodes.
 */

namespace kaori {

/**
 * @brief Handler for single-end combinatorial barcodes.
 *
 * In this design, the barcoding element is created from a template with multiple variable regions.
 * Each region contains a barcode from a different pool of options, where combinations are assembled randomly by library construction.
 * The idea is to use the large number of combinations to provide many unique identifiers, e.g., for cell-tracing applications.
 * This handler will capture the frequencies of each barcode combination. 
 *
 * @tparam max_size Maximum length of the template sequences on both reads.
 * @tparam num_variable Number of variable regions in the construct.
 */
template<size_t max_size, size_t num_variable>
class CombinatorialBarcodesSingleEnd {
public:
    /**
     * @brief Optional parameters for `CombinatorialBarcodeSingleEnd`.
     */
    struct Options {
        /**
         * Maximum number of mismatches allowed across the barcoding element.
         */
        int max_mismatches = 0;

        /** @param Whether to search only for the first match.
         * If `false`, the handler will search for the best match (i.e., fewest mismatches) instead.
         */
        bool use_first = true;

        /**
         * Strand(s) of the read sequence to search for the barcoding element.
         */
        SearchStrand strand = SearchStrand::FORWARD;

        /**
         * How duplicated barcode sequences should be handled.
         */
        DuplicateAction duplicates = DuplicateAction::ERROR;
    };

public:
    /**
     * @param[in] template_seq Template sequence for the first barcode.
     * This should contain exactly `num_variable` variable regions.
     * @param template_length Length of the template.
     * This should be less than or equal to `max_size`.
     * @param barcode_pools Array containing the known barcode sequences for each of the variable regions, in the order of their appearance in the template sequence.
     * @param options Optional parameters.
     *
     * @tparam BarcodePoolContainer Some iterable container of `BarcodePool` instances,
     * usually either a `std::vector` or a `std::array`.
     */
    template<class BarcodePoolContainer>
    CombinatorialBarcodesSingleEnd(const char* template_seq, size_t template_length, const BarcodePoolContainer& barcode_pools, const Options& options) :
        forward(search_forward(options.strand)),
        reverse(search_reverse(options.strand)),
        max_mm(options.max_mismatches),
        use_first(options.use_first),
        constant_matcher(template_seq, template_length, options.strand)
    {
        const auto& regions = constant_matcher.variable_regions();
        if (regions.size() != num_variable) { 
            throw std::runtime_error("expected " + std::to_string(num_variable) + " variable regions in the constant template");
        }
        if (barcode_pools.size() != num_variable) { 
            throw std::runtime_error("length of 'barcode_pools' should be equal to the number of variable regions");
        }

        for (size_t i = 0; i < num_variable; ++i) {
            size_t rlen = regions[i].second - regions[i].first;
            size_t vlen = barcode_pools[i].length;
            if (vlen != rlen) {
                throw std::runtime_error("length of variable region " + std::to_string(i + 1) + " (" + std::to_string(rlen) + 
                    ") should be the same as its sequences (" + std::to_string(vlen) + ")");
            }
        }

        // We'll be using this later.
        for (size_t i = 0; i < num_variable; ++i) {
            num_options[i] = barcode_pools[i].pool.size();
        }

        SimpleBarcodeSearch::Options bopt;
        bopt.max_mismatches = options.max_mismatches;
        bopt.duplicates = options.duplicates;

        if (forward) {
            bopt.reverse = false;
            for (size_t i = 0; i < num_variable; ++i) {
                forward_lib[i] = SimpleBarcodeSearch(barcode_pools[i], bopt);
            }
        }

        if (reverse) {
            bopt.reverse = true;
            for (size_t i = 0; i < num_variable; ++i) {
                reverse_lib[i] = SimpleBarcodeSearch(barcode_pools[num_variable - i - 1], bopt);
            }
        }
    }
        
    /**
     * @param t Whether to search only for the first match.
     * If `false`, the handler will search for the best match (i.e., fewest mismatches) instead.
     *
     * @return A reference to this `CombinatorialBarcodesSingleEnd` instance.
     */
    CombinatorialBarcodesSingleEnd& set_first(bool t = true) {
        use_first = t;
        return *this;
    }

public:
    /**
     * @cond
     */
    struct State {
        std::vector<std::array<int, num_variable> >collected;
        int total = 0;

        std::array<int, num_variable> temp;
        std::string buffer;

        // Default constructors should be called in this case, so it should be fine.
        std::array<typename SimpleBarcodeSearch::State, num_variable> forward_details, reverse_details;
    };
    /**
     * @endcond
     */

private:
    template<bool reverse>
    std::pair<bool, int> find_match(
        const char* seq, 
        size_t position, 
        int obs_mismatches, 
        const std::array<SimpleBarcodeSearch, num_variable>& libs, 
        std::array<typename SimpleBarcodeSearch::State, num_variable>& states, 
        std::array<int, num_variable>& temp,
        std::string& buffer
    ) const {
        const auto& regions = constant_matcher.template variable_regions<reverse>();

        for (size_t r = 0; r < num_variable; ++r) {
            auto range = regions[r];
            auto start = seq + position;
            buffer.clear(); // clear and insert preserves buffer's existing heap allocation.
            buffer.insert(buffer.end(), start + range.first, start + range.second);

            auto& curstate = states[r];
            libs[r].search(buffer, curstate, max_mm - obs_mismatches);
            if (curstate.index < 0) {
                return std::make_pair(false, 0);
            }
            
            obs_mismatches += curstate.mismatches;
            if (obs_mismatches > max_mm) {
                return std::make_pair(false, 0);
            }

            if constexpr(reverse) {
                temp[num_variable - r - 1] = curstate.index;
            } else {
                temp[r] = curstate.index;
            }
        }

        return std::make_pair(true, obs_mismatches);
    }

    std::pair<bool, int> forward_match(const char* seq, const typename ScanTemplate<max_size>::State& deets, State& state) const {
        return find_match<false>(seq, deets.position, deets.forward_mismatches, forward_lib, state.forward_details, state.temp, state.buffer);
    }

    std::pair<bool, int> reverse_match(const char* seq, const typename ScanTemplate<max_size>::State& deets, State& state) const {
        return find_match<true>(seq, deets.position, deets.reverse_mismatches, reverse_lib, state.reverse_details, state.temp, state.buffer);
    }

private:
    void process_first(State& state, const std::pair<const char*, const char*>& x) const {
        auto deets = constant_matcher.initialize(x.first, x.second - x.first);

        while (!deets.finished) {
            constant_matcher.next(deets);

            if (forward && deets.forward_mismatches <= max_mm) {
                if (forward_match(x.first, deets, state).first) {
                    state.collected.push_back(state.temp);
                    return;
                }
            }

            if (reverse && deets.reverse_mismatches <= max_mm) {
                if (reverse_match(x.first, deets, state).first) {
                    state.collected.push_back(state.temp);
                    return;
                }
            }
        }
    }

    void process_best(State& state, const std::pair<const char*, const char*>& x) const {
        auto deets = constant_matcher.initialize(x.first, x.second - x.first);
        bool found = false;
        int best_mismatches = max_mm + 1;
        std::array<int, num_variable> best_id;

        auto update = [&](std::pair<bool, int> match) -> void {
            if (match.first && match.second <= best_mismatches) {
                if (match.second == best_mismatches) {
                    if (best_id != state.temp) { // ambiguous.
                        found = false;
                    }
                } else { 
                    // A further optimization at this point would be to narrow
                    // max_mm to the current 'best_mismatches'. But
                    // this probably isn't worth it.

                    found = true;
                    best_mismatches = match.second;
                    best_id = state.temp;
                }
            }
        };

        while (!deets.finished) {
            constant_matcher.next(deets);

            if (forward && deets.forward_mismatches <= max_mm) {
                update(forward_match(x.first, deets, state));
            }

            if (reverse && deets.reverse_mismatches <= max_mm) {
                update(reverse_match(x.first, deets, state));
            }
        }

        if (found) {
            state.collected.push_back(best_id);
        }
    }

public:
    /**
     * @cond
     */
    State initialize() const {
        return State();
    }

    void reduce(State& s) {
        if (forward) {
            for (size_t r = 0; r < num_variable; ++r) {
                forward_lib[r].reduce(s.forward_details[r]);
            }
        }
        if (reverse) {
            for (size_t r = 0; r < num_variable; ++r) {
                reverse_lib[r].reduce(s.reverse_details[r]);
            }
        }

        combinations.insert(combinations.end(), s.collected.begin(), s.collected.end());
        total += s.total;
        return;
    }

    void process(State& state, const std::pair<const char*, const char*>& x) const {
        if (use_first) {
            process_first(state, x);
        } else {
            process_best(state, x);
        }
        ++state.total;
    }

    static constexpr bool use_names = false;
    /**
     * @endcond
     */

public:
    /**
     * Sort the combinations for easier frequency counting.
     */
    void sort() {
        sort_combinations(combinations, num_options);
    }

    /**
     * @return All combinations encountered by the handler.
     */
    const std::vector<std::array<int, num_variable> >& get_combinations() const {
        return combinations;
    }

    /**
     * @return Total number of reads processed by the handler.
     */
    int get_total() const {
        return total;
    }
private:
    bool forward;
    bool reverse;
    int max_mm;
    bool use_first;
    size_t nregions;

    ScanTemplate<max_size> constant_matcher;
    std::array<SimpleBarcodeSearch, num_variable> forward_lib, reverse_lib;
    std::array<size_t, num_variable> num_options;

    std::vector<std::array<int, num_variable> > combinations;
    int total = 0;
};

}

#endif
