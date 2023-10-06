#ifndef KAORI_DUAL_BARCODES_SINGLE_END_HPP
#define KAORI_DUAL_BARCODES_SINGLE_END_HPP

#include "../ScanTemplate.hpp"
#include "../BarcodeSearch.hpp"
#include "../utils.hpp"

#include <array>
#include <vector>

/**
 * @file DualBarcodesSingleEnd.hpp
 *
 * @brief Process single-end dual barcodes.
 */

namespace kaori {

/**
 * @brief Handler for single-end dual barcodes.
 *
 * In this design, the barcoding element is created from a template with multiple variable regions.
 * Each region contains a barcode from a different pool of options, where the valid combinations of barcodes across variable regions are known beforehand.
 * This differs from `CombinatorialBarcodesSingleEnd` where the combinations are assembled randomly.
 * Despite its name, this handler can actually handle any number (>= 2) of variable regions in the combination.
 * It will count the frequency of each barcode combination, along with the total number of reads. 
 *
 * @tparam max_size Maximum length of the template sequence.
 */
template<size_t max_size>
class DualBarcodesSingleEnd {
public:
    /**
     * @brief Optional parameters for `DualBarcodeSingleEnd`.
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
     * @param[in] template_seq Template sequence containing any number (usually 2 or more) of variable regions.
     * @param template_length Length of the template.
     * This should be less than or equal to `max_size`.
     * @param barcode_pools Array containing the known barcode sequences for each of the variable regions, in the order of their appearance in the template sequence.
     * Each pool should have the same length, and corresponding values across pools define a specific combination of barcodes. 
     * @param options Optional parameters.
     */
    DualBarcodesSingleEnd(const char* template_seq, size_t template_length, const std::vector<BarcodePool>& barcode_pools, const Options& options) :
        forward(search_forward(options.strand)),
        reverse(search_reverse(options.strand)),
        max_mm(options.max_mismatches),
        use_first(options.use_first),
        constant_matcher(template_seq, template_length, options.strand)
    {
        const auto& regions = constant_matcher.variable_regions();
        num_variable = regions.size();
        if (barcode_pools.size() != num_variable) {
            throw std::runtime_error("length of 'barcode_pools' should equal the number of variable regions");
        }

        for (size_t i = 0; i < num_variable; ++i) {
            size_t rlen = regions[i].second - regions[i].first;
            size_t vlen = barcode_pools[i].length;
            if (vlen != rlen) {
                throw std::runtime_error("length of variable region " + std::to_string(i + 1) + " (" + std::to_string(rlen) + 
                    ") should be the same as its sequences (" + std::to_string(vlen) + ")");
            }
        }

        size_t num_choices = 0;
        if (num_variable) {
            num_choices = barcode_pools[0].size();
            for (size_t i = 1; i < num_variable; ++i) {
                if (num_choices != barcode_pools[i].size()) {
                    throw std::runtime_error("all entries of 'barcode_pools' should have the same length");
                }
            }
        }
        counts.resize(num_choices);

        // Constructing the combined varlib.
        std::vector<std::string> combined(num_choices); 
        for (size_t v = 0; v < num_variable; ++v) {
            const auto& curpool = barcode_pools[v];
            size_t n = curpool.length;
            for (size_t c = 0; c < num_choices; ++c) {
                auto ptr = curpool[c];
                combined[c].insert(combined[c].end(), ptr, ptr + n);
            }
        }

        SimpleBarcodeSearch::Options bopt;
        bopt.max_mismatches = options.max_mismatches;
        bopt.duplicates = options.duplicates;

        if (forward) {
            bopt.reverse = false;
            forward_lib = SimpleBarcodeSearch(combined, bopt);
        }

        if (reverse) {
            bopt.reverse = true;
            reverse_lib = SimpleBarcodeSearch(combined, bopt);
        }
    }

public:
    /**
     * @cond
     */
    struct State {
        std::vector<int> counts;
        int total = 0;

        std::string buffer;

        // Default constructors should be called in this case, so it should be fine.
        typename SimpleBarcodeSearch::State forward_details, reverse_details;
    };
    /**
     * @endcond
     */

private:
    template<bool reverse>
    std::pair<int, int> find_match(
        const char* seq, 
        size_t position, 
        int obs_mismatches, 
        const SimpleBarcodeSearch& lib, 
        typename SimpleBarcodeSearch::State state, 
        std::string& buffer
    ) const {
        const auto& regions = constant_matcher.template variable_regions<reverse>();
        buffer.clear();

        for (size_t r = 0; r < num_variable; ++r) {
            auto start = seq + position;
            buffer.insert(buffer.end(), start + regions[r].first, start + regions[r].second);
        }

        lib.search(buffer, state, max_mm - obs_mismatches);
        return std::make_pair(state.index, obs_mismatches + state.mismatches);
    }

    std::pair<int, int> forward_match(const char* seq, const typename ScanTemplate<max_size>::State& deets, State& state) const {
        return find_match<false>(seq, deets.position, deets.forward_mismatches, forward_lib, state.forward_details, state.buffer);
    }

    std::pair<int, int> reverse_match(const char* seq, const typename ScanTemplate<max_size>::State& deets, State& state) const {
        return find_match<true>(seq, deets.position, deets.reverse_mismatches, reverse_lib, state.reverse_details, state.buffer);
    }

private:
    bool process_first(State& state, const std::pair<const char*, const char*>& x) const {
        auto deets = constant_matcher.initialize(x.first, x.second - x.first);

        while (!deets.finished) {
            constant_matcher.next(deets);

            if (forward && deets.forward_mismatches <= max_mm) {
                auto id = forward_match(x.first, deets, state).first;
                if (id >= 0) {
                    ++state.counts[id];
                    return true;
                }
            }

            if (reverse && deets.reverse_mismatches <= max_mm) {
                auto id = reverse_match(x.first, deets, state).first;
                if (id >= 0) {
                    ++state.counts[id];
                    return true;
                }
            }
        }
        return false;
    }

    bool process_best(State& state, const std::pair<const char*, const char*>& x) const {
        auto deets = constant_matcher.initialize(x.first, x.second - x.first);
        bool found = false;
        int best_mismatches = max_mm + 1;
        int best_id = -1;

        auto update = [&](std::pair<int, int> match) -> void {
            if (match.first < 0){ 
                return;
            }
            if (match.second == best_mismatches) {
                if (best_id != match.first) { // ambiguous.
                    found = false;
                }
            } else if (match.second < best_mismatches) { 
                // A further optimization at this point would be to narrow
                // max_mm to the current 'best_mismatches'. But
                // this probably isn't worth it.

                found = true;
                best_mismatches = match.second;
                best_id = match.first;
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
            ++state.counts[best_id];
        }
        return found;
    }

public:
    /**
     * @cond
     */
    State initialize() const {
        State output;
        output.counts.resize(counts.size());
        return output;
    }

    void reduce(State& s) {
        if (forward) {
            forward_lib.reduce(s.forward_details);
        }
        if (reverse) {
            reverse_lib.reduce(s.reverse_details);
        }

        for (size_t i = 0, end = counts.size(); i < end; ++i) {
            counts[i] += s.counts[i];
        }
        total += s.total;
        return;
    }

    bool process(State& state, const std::pair<const char*, const char*>& x) const {
        ++state.total;
        if (use_first) {
            return process_first(state, x);
        } else {
            return process_best(state, x);
        }
    }

    static constexpr bool use_names = false;
    /**
     * @endcond
     */

public:
    /**
     * @return Counts for each combination.
     */
    const std::vector<int>& get_counts() const {
        return counts;
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

    ScanTemplate<max_size> constant_matcher;
    size_t num_variable;

    SimpleBarcodeSearch forward_lib, reverse_lib;
    std::vector<int> counts;
    int total = 0;
};

}

#endif
