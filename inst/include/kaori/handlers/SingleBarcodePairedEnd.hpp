#ifndef KAORI_SINGLE_BARCODE_PAIRED_END_HPP
#define KAORI_SINGLE_BARCODE_PAIRED_END_HPP

#include "../SimpleSingleMatch.hpp"
#include <vector>

/**
 * @file SingleBarcodePairedEnd.hpp
 *
 * @brief Process paired-end single barcodes.
 */

namespace kaori {

/**
 * @brief Handler for paired-end single barcodes.
 *
 * In this design, the barcoding element is created from a template with a single variable region drawn from a pool of barcode sequences.
 * The construct containing the barcoding element is then subjected to paired-end sequencing, where either end could contain the barcoding element.
 * This handler will search both reads for the barcoding element and count the frequency of each barcode.
 *
 * @tparam max_size Maximum length of the template sequence.
 */
template<size_t max_size>
class SingleBarcodePairedEnd {
public:
    /**
     * @brief Optional parameters for `SingleBarcodePairedEnd`.
     */
    struct Options {
        /** 
         * Whether to search only for the first match.
         * If `false`, the handler will search for the best match (i.e., fewest mismatches) instead.
         */
        bool use_first = true;

        /** 
         * Maximum number of mismatches allowed across the barcoding element.
         */
        int max_mismatches = 0;

        /**
         * Strand(s) of each read sequence to search.
         */
        SearchStrand strand = SearchStrand::FORWARD;

        /** 
         * How duplicated barcode sequences should be handled.
         */
        DuplicateAction duplicates = DuplicateAction::ERROR;
    };

public:
    /**
     * @param[in] template_seq Template sequence of the barcoding element.
     * This should contain exactly one variable region.
     * @param template_length Length of the template.
     * This should be less than or equal to `max_size`.
     * @param barcode_pool Known barcode sequences for the variable region.
     * @param options Optional parameters.
     */
    SingleBarcodePairedEnd(const char* template_seq, size_t template_length, const BarcodePool& barcode_pool, const Options& options) :
        matcher(
            template_seq, 
            template_length, 
            barcode_pool, 
            [&]{
                typename SimpleSingleMatch<max_size>::Options ssopt;
                ssopt.strand = options.strand;
                ssopt.max_mismatches = options.max_mismatches;
                ssopt.duplicates = options.duplicates;
                return ssopt;
            }()
        ),
        counts(barcode_pool.size()),
        use_first(options.use_first)
    {}

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

    void process(State& state, const std::pair<const char*, const char*>& r1, const std::pair<const char*, const char*>& r2) const {
        if (use_first) {
            if (matcher.search_first(r1.first, r1.second - r1.first, state.search)) {
                ++state.counts[state.search.index];
            } else if (matcher.search_first(r2.first, r2.second - r2.first, state.search)) {
                ++state.counts[state.search.index];
            }
        } else {
            bool found1 = matcher.search_best(r1.first, r1.second - r1.first, state.search);
            auto id1 = state.search.index;
            auto mm1 = state.search.mismatches;

            bool found2 = matcher.search_best(r2.first, r2.second - r2.first, state.search);
            auto id2 = state.search.index;
            auto mm2 = state.search.mismatches;

            if (found1 && !found2) {
                ++state.counts[id1];
            } else if (!found1 && found2) {
                ++state.counts[id2];
            } else if (found1 && found2) {
                if (mm1 < mm2) {
                    ++state.counts[id1];
                } else if (mm1 > mm2) {
                    ++state.counts[id2];
                } else if (id1 == id2) {
                    ++state.counts[id1];
                }
            }
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
