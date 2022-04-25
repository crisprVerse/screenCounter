#ifndef KAORI_DUAL_BARCODES_WITH_DIAGNOSTICS_HPP
#define KAORI_DUAL_BARCODES_WITH_DIAGNOSTICS_HPP

#include "DualBarcodes.hpp"
#include "CombinatorialBarcodesPairedEnd.hpp"
#include "../utils.hpp"

/**
 * @file DualBarcodesWithDiagnostics.hpp
 *
 * @brief Process dual barcodes with extra diagnostics.
 */

namespace kaori {

/**
 * @brief Handler for dual barcodes with extra diagnostics.
 *
 * This provides the same information as `DualBarcodes` but also captures the frequency of the invalid combinations.
 * These frequences can be helpful for diagnosing problems with library construction.
 * The handler also counts the number of reads where only one barcode construct matches to a read.
 *
 * @tparam max_size Maximum length of the template sequences on both reads.
 */
template<size_t max_size>
class DualBarcodesWithDiagnostics { 
public:
    /**
     * @param[in] template_seq1 Pointer to a character array containing the first template sequence. 
     * This should contain exactly one variable region.
     * @param template_length1 Length of the first template.
     * This should be less than or equal to `max_size`.
     * @param reverse1 Whether to search the reverse strand of the read for the first template.
     * @param barcode_pool1 Pool of known barcode sequences for the variable region in the first template.
     * @param max_mismatches1 Maximum number of mismatches across the target sequence corresponding to the first template.
     * @param[in] template_seq2 Pointer to a character array containing the second template sequence. 
     * This should contain exactly one variable region.
     * @param template_length2 Length of the second template.
     * This should be less than or equal to `max_size`.
     * @param reverse2 Whether to search the reverse strand of the read for the second template.
     * @param barcode_pool2 Pool of known barcode sequences for the variable region in the second template.
     * @param max_mismatches2 Maximum number of mismatches across the target sequence corresponding to the second template.
     * @param random Whether the reads are randomized with respect to the first/second target sequences.
     * If `false`, the first read is searched for the first target sequence only, and the second read is searched for the second target sequence only.
     * If `true`, an additional search will be performed in the opposite orientation.
     *
     * `barcode_pool1` and `barcode_pool2` are expected to have the same number of barcodes (possibly duplicated).
     * Corresponding values across the two pools define a particular combination of dual barcodes. 
     */
    DualBarcodesWithDiagnostics(
        const char* template_seq1, size_t template_length1, bool reverse1, const BarcodePool& barcode_pool1, int max_mismatches1, 
        const char* template_seq2, size_t template_length2, bool reverse2, const BarcodePool& barcode_pool2, int max_mismatches2,
        bool random = false
    ) :
        dual_handler(template_seq1, template_length1, reverse1, barcode_pool1, max_mismatches1, template_seq2, template_length2, reverse2, barcode_pool2, max_mismatches2, random),

        // we allow duplicates in the trie.
        combo_handler(template_seq1, template_length1, reverse1, barcode_pool1, max_mismatches1, template_seq2, template_length2, reverse2, barcode_pool2, max_mismatches2, random, true) 
    {}

    /**
     * @param t Whether to search only for the first match across reads (for valid combinations) or in each read (for invalid combinations).
     * If `false`, the handler will search for the best match (i.e., fewest mismatches) instead.
     *
     * @return A reference to this `DualBarcodesWithDiagnostics` instance.
     */
    DualBarcodesWithDiagnostics& set_first(bool t = true) {
        dual_handler.set_first(t);
        combo_handler.set_first(t);
        return *this;
    }

private:
    DualBarcodes<max_size> dual_handler;
    CombinatorialBarcodesPairedEnd<max_size> combo_handler;

public:
    /**
     *@cond
     */
    struct State {
        State() {}
        State(typename DualBarcodes<max_size>::State ds, typename CombinatorialBarcodesPairedEnd<max_size>::State cs) : dual_state(std::move(ds)), combo_state(std::move(cs)) {}

        /**
         * @cond
         */
        typename DualBarcodes<max_size>::State dual_state;
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
     * @return Sort the invalid combinations for easier frequency counting.
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

}

#endif
