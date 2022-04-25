#ifndef KAORI_SINGLE_BARCODE_SINGLE_END_HPP
#define KAORI_SINGLE_BARCODE_SINGLE_END_HPP

#include "../SimpleSingleMatch.hpp"
#include <vector>

/**
 * @file SingleBarcodeSingleEnd.hpp
 *
 * @brief Process single-end single barcodes.
 */

namespace kaori {

/**
 * @brief Handler for single-end single barcodes.
 *
 * In this design, the target sequence is created from a template with a single variable region drawn from a pool of barcode sequences.
 * The construct containing the target sequence is then subjected to single-end sequencing.
 * This handler will search the read for the target sequence and count the frequency of each barcode.
 *
 * @tparam max_size Maximum length of the template sequences on both reads.
 */
template<size_t max_size>
class SingleBarcodeSingleEnd {
public:
    /**
     * @param[in] template_seq Template sequence for the first barcode.
     * This should contain exactly one variable region.
     * @param template_length Length of the template.
     * This should be less than or equal to `max_size`.
     * @param strand Strand to use for searching the read sequence - forward (0), reverse (1) or both (2).
     * @param barcode_pool Known barcode sequences for the variable region.
     * @param max_mismatches Maximum number of mismatches allowed across the target sequence.
     */
    SingleBarcodeSingleEnd(const char* template_seq, size_t template_length, int strand, const BarcodePool& barcode_pool, int max_mismatches = 0) : 
        matcher(template_seq, template_length, strand != 1, strand != 0, barcode_pool, max_mismatches), counts(barcode_pool.size()) {}

    /**
     * @param t Whether to search only for the first match.
     * If `false`, the handler will search for the best match (i.e., fewest mismatches) instead.
     *
     * @return A reference to this `SingleBarcodeSingleEnd` instance.
     */
    SingleBarcodeSingleEnd& set_first(bool t = true) {
        use_first = t;
        return *this;
    }

public:
    /**
     * @cond
     */
    struct State {
        State() {}

        State(typename SimpleSingleMatch<max_size>::State s, size_t nvar) : search(std::move(s)), counts(nvar) {}

        typename SimpleSingleMatch<max_size>::State search;
        std::vector<int> counts;
        int total = 0;
    };

    void process(State& state, const std::pair<const char*, const char*>& x) const {
        bool found = false;
        if (use_first) {
            found = matcher.search_first(x.first, x.second - x.first, state.search);
        } else {
            found = matcher.search_best(x.first, x.second - x.first, state.search);
        }
        if (found) {
            ++(state.counts[state.search.index]);
        }
        ++state.total;
    }

    static constexpr bool use_names = false;
    /**
     * @endcond
     */

public:
    /**
     * @cond
     */
    State initialize() const {
        return State(matcher.initialize(), counts.size());
    }

    void reduce(State& s) {
        matcher.reduce(s.search);
        for (size_t i = 0; i < counts.size(); ++i) {
            counts[i] += s.counts[i];
        }
        total += s.total;
    }
    /**
     * @endcond
     */

private:
    SimpleSingleMatch<max_size> matcher;
    std::vector<int> counts;
    int total = 0;
    bool use_first = true;

public:
    /**
     * @return Vector containing the frequency of each barcode.
     * This has length equal to the number of valid barcodes (i.e., the length of `barcode_pool` in the constructor).
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
};

}

#endif
