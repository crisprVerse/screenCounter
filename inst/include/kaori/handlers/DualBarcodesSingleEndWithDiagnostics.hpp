#ifndef KAORI_DUAL_BARCODES_SINGLE_END_WITH_DIAGNOSTICS_HPP
#define KAORI_DUAL_BARCODES_SINGLE_END_WITH_DIAGNOSTICS_HPP

#include "DualBarcodesSingleEnd.hpp"
#include "CombinatorialBarcodesSingleEnd.hpp"
#include "../utils.hpp"

/**
 * @file DualBarcodesSingleEndWithDiagnostics.hpp
 *
 * @brief Process dual barcodes with extra diagnostics.
 */

namespace kaori {

/**
 * @brief Handler for dual barcodes with extra diagnostics.
 *
 * This provides the same information as `DualBarcodesSingleEnd` but also captures the frequency of the invalid combinations.
 * These frequences can be helpful for diagnosing problems with library construction.
 *
 * @tparam max_size Maximum length of the template sequences on both reads.
 * @tparam num_variable Number of the template sequences on both reads.
 */
template<size_t max_size, size_t num_variable>
class DualBarcodesSingleEndWithDiagnostics { 
public:
    /**
     * @param[in] template_seq Pointer to a character array containing the template sequence. 
     * @param template_length Length of the template.
     * This should be less than or equal to `max_size`.
     * @param barcode_pools Pools of known barcode sequences for each variable region in the template.
     * Each pool should have the same length, and corresponding values across pools define a specific combination of barcodes. 
     * @param options Optional parameters.
     */
    DualBarcodesSingleEndWithDiagnostics(
        const char* template_seq, 
        size_t template_length, 
        const std::vector<BarcodePool>& barcode_pools, 
        const typename DualBarcodesSingleEnd<max_size>::Options& options
    ) :
        dual_handler(template_seq, template_length, barcode_pools, options),

        combo_handler(template_seq, template_length, barcode_pools, 
            [&]{
                typename CombinatorialBarcodesSingleEnd<max_size, num_variable>::Options combopt;
                combopt.use_first = options.use_first;

                combopt.max_mismatches = options.max_mismatches;
                combopt.strand = options.strand;

                // we allow duplicates in the trie for each individual barcode, as only the combinations are unique in the dual barcode setup.
                combopt.duplicates = DuplicateAction::FIRST; 
                return combopt;
            }()
        )
    {}

private:
    DualBarcodesSingleEnd<max_size> dual_handler;
    CombinatorialBarcodesSingleEnd<max_size, num_variable> combo_handler;

public:
    /**
     *@cond
     */
    struct State {
        State() {}
        State(typename DualBarcodesSingleEnd<max_size>::State ds, typename CombinatorialBarcodesSingleEnd<max_size, num_variable>::State cs) : dual_state(std::move(ds)), combo_state(std::move(cs)) {}

        /**
         * @cond
         */
        typename DualBarcodesSingleEnd<max_size>::State dual_state;
        typename CombinatorialBarcodesSingleEnd<max_size, num_variable>::State combo_state;
        /**
         * @endcond
         */
    };

    State initialize() const {
        return State(dual_handler.initialize(), combo_handler.initialize());
    }

    void reduce(State& s) {
        dual_handler.reduce(s.dual_state);
        combo_handler.reduce(s.combo_state);
    }

    constexpr static bool use_names = false;
    /**
     *@endcond
     */

public:
    /**
     *@cond
     */
    void process(State& state, const std::pair<const char*, const char*>& x) const {
        // Only searching for combinations if we couldn't find a proper dual barcode match.
        if (!dual_handler.process(state.dual_state, x)) {
            combo_handler.process(state.combo_state, x);
        }
    }
    /**
     *@endcond
     */

public:
    /**
     * Sort the invalid combinations for easier frequency counting.
     * Combinations are sorted by the first index, and then the second index.
     */
    void sort() {
        combo_handler.sort();
    }

    /**
     * @return Vector containing the frequency of each valid combination.
     * This has length equal to the number of valid dual barcode combinations (i.e., the length of `barcode_pool1` and `barcode_pool2` in the constructor).
     * Each entry contains the count for the corresponding dual barcode combination.
     */
    const std::vector<int>& get_counts() const {
        return dual_handler.get_counts();
    }

    /**
     * @return All invalid combinations encountered by the handler.
     * In each array, the first and second element contains the indices of known barcodes in the first and second pools, respectively.
     */
    const std::vector<std::array<int, num_variable> >& get_combinations() const {
        return combo_handler.get_combinations();
    }

    /**
     * @return Total number of reads processed by the handler.
     */
    int get_total() const {
        return dual_handler.get_total();
    }
};

}

#endif
