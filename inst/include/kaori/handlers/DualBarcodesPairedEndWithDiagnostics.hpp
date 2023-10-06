#ifndef KAORI_DUAL_BARCODES_PAIRED_END_WITH_DIAGNOSTICS_HPP
#define KAORI_DUAL_BARCODES_PAIRED_END_WITH_DIAGNOSTICS_HPP

#include "DualBarcodesPairedEnd.hpp"
#include "CombinatorialBarcodesPairedEnd.hpp"
#include "../utils.hpp"

/**
 * @file DualBarcodesPairedEndWithDiagnostics.hpp
 *
 * @brief Process dual barcodes with extra diagnostics.
 */

namespace kaori {

/**
 * @brief Handler for dual barcodes with extra diagnostics.
 *
 * This provides the same information as `DualBarcodesPairedEnd` but also captures the frequency of the invalid combinations.
 * These frequences can be helpful for diagnosing problems with library construction.
 * The handler also counts the number of reads where only one barcode construct matches to a read.
 *
 * @tparam max_size Maximum length of the template sequences on both reads.
 */
template<size_t max_size>
class DualBarcodesPairedEndWithDiagnostics { 
public:
    /**
     * @param[in] template_seq1 Pointer to a character array containing the first template sequence. 
     * This should contain exactly one variable region.
     * @param template_length1 Length of the first template.
     * This should be less than or equal to `max_size`.
     * @param barcode_pool1 Pool of known barcode sequences for the variable region in the first template.
     * @param[in] template_seq2 Pointer to a character array containing the second template sequence. 
     * This should contain exactly one variable region.
     * @param template_length2 Length of the second template.
     * This should be less than or equal to `max_size`.
     * @param barcode_pool2 Pool of known barcode sequences for the variable region in the second template.
     * @param options Optional parameters.
     *
     * `barcode_pool1` and `barcode_pool2` are expected to have the same number of barcodes (possibly duplicated).
     * Corresponding values across the two pools define a particular combination of dual barcodes. 
     */
    DualBarcodesPairedEndWithDiagnostics(
        const char* template_seq1, size_t template_length1, const BarcodePool& barcode_pool1,
        const char* template_seq2, size_t template_length2, const BarcodePool& barcode_pool2, 
        const typename DualBarcodesPairedEnd<max_size>::Options& options
    ) :
        dual_handler(template_seq1, template_length1, barcode_pool1, template_seq2, template_length2, barcode_pool2, options),

        combo_handler(
            template_seq1, 
            template_length1, 
            barcode_pool1, 
            template_seq2, 
            template_length2, 
            barcode_pool2, 
            [&]{
                typename CombinatorialBarcodesPairedEnd<max_size>::Options combopt;
                combopt.use_first = options.use_first;

                combopt.max_mismatches1 = options.max_mismatches1;
                combopt.strand1 = options.strand1;
                combopt.max_mismatches2 = options.max_mismatches2;
                combopt.strand2 = options.strand2;

                // we allow duplicates in the trie for each individual barcode, as only the pairs are unique in the dual barcode setup.
                combopt.duplicates = DuplicateAction::FIRST; 
                combopt.random = options.random;
                return combopt;
            }()
        )
    {}

private:
    DualBarcodesPairedEnd<max_size> dual_handler;
    CombinatorialBarcodesPairedEnd<max_size> combo_handler;

public:
    /**
     *@cond
     */
    struct State {
        State() {}
        State(typename DualBarcodesPairedEnd<max_size>::State ds, typename CombinatorialBarcodesPairedEnd<max_size>::State cs) : dual_state(std::move(ds)), combo_state(std::move(cs)) {}

        /**
         * @cond
         */
        typename DualBarcodesPairedEnd<max_size>::State dual_state;
        typename CombinatorialBarcodesPairedEnd<max_size>::State combo_state;
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
    void process(State& state, const std::pair<const char*, const char*>& r1, const std::pair<const char*, const char*>& r2) const {
        // Only searching for combinations if we couldn't find a proper dual barcode match.
        if (!dual_handler.process(state.dual_state, r1, r2)) {
            combo_handler.process(state.combo_state, r1, r2);
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
    const std::vector<std::array<int, 2> >& get_combinations() const {
        return combo_handler.get_combinations();
    }

    /**
     * @return Total number of read pairs processed by the handler.
     */
    int get_total() const {
        return dual_handler.get_total();
    }

    /**
     * @return Number of read pairs with a valid match to the first barcode but no valid match to the second barcode.
     */
    int get_barcode1_only() const {
        return combo_handler.get_barcode1_only();
    }

    /**
     * @return Number of read pairs with a valid match to the second barcode but no valid match to the first barcode.
     */
    int get_barcode2_only() const {
        return combo_handler.get_barcode2_only();
    }
};

/**
 * @cond
 */
// Soft-deprecated back-compatible aliases.
template<size_t max_size>
using DualBarcodesWithDiagnostics = DualBarcodesPairedEndWithDiagnostics<max_size>;
/**
 * @endcond
 */

}

#endif
