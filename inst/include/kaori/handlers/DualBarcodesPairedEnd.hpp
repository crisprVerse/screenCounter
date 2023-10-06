#ifndef KAORI_DUAL_BARCODES_HPP
#define KAORI_DUAL_BARCODES_HPP

#include "../ScanTemplate.hpp"
#include "../BarcodeSearch.hpp"
#include "../utils.hpp"

/**
 * @file DualBarcodesPairedEnd.hpp
 *
 * @brief Process dual barcodes.
 */

namespace kaori {

/**
 * @brief Handler for dual barcodes.
 *
 * In this design, each read contains a barcoding element created from a template with a single variable region.
 * For one read, the barcode is drawn from one pool of options, while the other read contains a barcode from another pool.
 * However, unlike `CombinatorialBarcodesPairedEnd`, the combinations are not random but are specifically assembled, typically corresponding to specific pairs of genes.
 * This handler will capture the frequencies of each barcode combination. 
 *
 * @tparam max_size Maximum length of the template sequences on both reads.
 */
template<size_t max_size>
class DualBarcodesPairedEnd { 
public:
    /**
     * @brief Optional parameters for `DualBarcodesPairedEnd`.
     */
    struct Options {
        /** 
         * Whether to search only for the first match.
         * If `false`, the handler will search for the best match (i.e., fewest mismatches) instead.
         */
        bool use_first = true;

        /** 
         * Maximum number of mismatches allowed across the first barcoding element.
         */
        int max_mismatches1 = 0;

        /**
         * Strand of the read sequence to search for the first barcoding element.
         * `BOTH` is not supported right now... sorry.
         */
        SearchStrand strand1 = SearchStrand::FORWARD;

        /** 
         * Maximum number of mismatches allowed across the second barcoding element.
         */
        int max_mismatches2 = 0;

        /**
         * Strand of the read sequence to search for the second barcoding element.
         * `BOTH` is not supported right now... sorry.
         */
        SearchStrand strand2 = SearchStrand::FORWARD;

        /**
         * How duplicated pairs of barcode sequences should be handled.
         */
        DuplicateAction duplicates = DuplicateAction::ERROR;

        /**
         * Whether the reads are randomized with respect to the first/second barcoding elements.
         * If `false`, the first read is searched for the first barcoding element only, and the second read is searched for the second barcoding element only.
         * If `true`, an additional search will be performed in the opposite orientation.
         */
        bool random = false;
    };

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
     * `barcode_pool1` and `barcode_pool2` are expected to have the same number of barcodes.
     * Corresponding values across the two pools define a particular combination of dual barcodes. 
     * Duplication of sequences within each pool is allowed; only pairs of the same barcodes are considered to be duplicates with respect to `Options::duplicates`.
     */
    DualBarcodesPairedEnd(
        const char* template_seq1, size_t template_length1, const BarcodePool& barcode_pool1, 
        const char* template_seq2, size_t template_length2, const BarcodePool& barcode_pool2,
        const Options& options
    ) :
        search_reverse1(search_reverse(options.strand1)),
        search_reverse2(search_reverse(options.strand2)),
        constant1(template_seq1, template_length1, options.strand1),
        constant2(template_seq2, template_length2, options.strand2),
        max_mm1(options.max_mismatches1),
        max_mm2(options.max_mismatches2),
        randomized(options.random),
        use_first(options.use_first)
    {
        auto num_options = barcode_pool1.size();
        if (num_options != barcode_pool2.size()) {
            throw std::runtime_error("both barcode pools should be of the same length");
        }
        counts.resize(num_options);

        size_t len1;
        {
            const auto& regions = constant1.variable_regions();
            if (regions.size() != 1) { 
                throw std::runtime_error("expected one variable region in the first constant template");
            }
            len1 = regions[0].second - regions[0].first;
            if (len1 != barcode_pool1.length) {
                throw std::runtime_error("length of variable sequences (" + std::to_string(barcode_pool1.length) + 
                    ") should be the same as the variable region (" + std::to_string(len1) + ")");
            }
        }

        size_t len2;
        {
            const auto& regions = constant2.variable_regions();
            if (regions.size() != 1) { 
                throw std::runtime_error("expected one variable region in the second constant template");
            }
            len2 = regions[0].second - regions[0].first;
            if (len2 != barcode_pool2.length) {
                throw std::runtime_error("length of variable sequences (" + std::to_string(barcode_pool2.length) + 
                    ") should be the same as the variable region (" + std::to_string(len2) + ")");
            }
        }

        // Constructing the combined strings.
        std::vector<std::string> combined;
        combined.reserve(num_options);

        for (size_t i = 0; i < num_options; ++i) {
            std::string current;

            auto ptr1 = barcode_pool1[i];
            if (search_reverse1) {
                for (size_t j = 0; j < len1; ++j) {
                    current += complement_base<true, true>(ptr1[len1 - j - 1]);
                }
            } else {
                current.insert(current.end(), ptr1, ptr1 + len1);
            }

            auto ptr2 = barcode_pool2[i];
            if (search_reverse2) {
                for (size_t j = 0; j < len2; ++j) {
                    current += complement_base<true, true>(ptr2[len2 - j - 1]);
                }
            } else {
                current.insert(current.end(), ptr2, ptr2 + len2);
            }

            combined.push_back(std::move(current));
        }

        // Constructing the combined varlib.
        BarcodePool combined_set(combined);
        varlib = SegmentedBarcodeSearch<2>(
            combined_set,
            std::array<int, 2>{ static_cast<int>(len1), static_cast<int>(len2) }, 
            [&]{
                typename SegmentedBarcodeSearch<2>::Options bopt;
                bopt.max_mismatches = { max_mm1, max_mm2 };
                bopt.duplicates = options.duplicates;
                bopt.reverse = false; // we already handle strandedness when creating 'combined_set'.
                return bopt;
            }()
        );
    }

public:
    /**
     *@cond
     */
    struct State {
        State(size_t n = 0) : counts(n) {}
        std::vector<int> counts;
        int total = 0;

        std::pair<std::string, int> first_match;
        std::vector<std::pair<std::string, int> > second_matches;
        std::string combined;

        // Default constructors should be called in this case, so it should be fine.
        typename SegmentedBarcodeSearch<2>::State details;
    };

    State initialize() const {
        return State(counts.size());
    }

    void reduce(State& s) {
        varlib.reduce(s.details);
        for (size_t i = 0; i < counts.size(); ++i) {
            counts[i] += s.counts[i];
        }
        total += s.total;
    }

    constexpr static bool use_names = false;
    /**
     * @endcond
     */

private:
    static void fill_store(std::pair<std::string, int>& first_match, const char* start, const char* end, int mm) {
        first_match.first.clear();
        first_match.first.insert(first_match.first.end(), start, end);
        first_match.second = mm;
        return;
    }

    static void fill_store(std::vector<std::pair<std::string, int> >& second_matches, const char* start, const char* end, int mm) {
        second_matches.emplace_back(std::string(start, end), mm);
        return;
    }

    template<class Store>
    static bool inner_process(
        bool reverse, 
        const ScanTemplate<max_size>& constant, 
        int max_mm,
        const char* against,
        typename ScanTemplate<max_size>::State& deets,
        Store& store)
    {
        while (!deets.finished) {
            constant.next(deets);
            if (reverse) {
                if (deets.reverse_mismatches <= max_mm) {
                    const auto& reg = constant.template variable_regions<true>()[0];
                    auto start = against + deets.position;
                    fill_store(store, start + reg.first, start + reg.second, deets.reverse_mismatches);
                    return true;
                }
            } else {
                if (deets.forward_mismatches <= max_mm) {
                    const auto& reg = constant.variable_regions()[0];
                    auto start = against + deets.position;
                    fill_store(store, start + reg.first, start + reg.second, deets.forward_mismatches);
                    return true;
                }
            }
        }
        return false;
    }

    bool process_first(State& state, const std::pair<const char*, const char*>& against1, const std::pair<const char*, const char*>& against2) const {
        auto deets1 = constant1.initialize(against1.first, against1.second - against1.first);
        auto deets2 = constant2.initialize(against2.first, against2.second - against2.first);

        state.second_matches.clear();

        auto checker = [&](size_t idx2) -> bool {
            const auto& current2 = state.second_matches[idx2];
            state.combined = state.first_match.first;
            state.combined += current2.first; // on a separate line to avoid creating a std::string intermediate.
            varlib.search(state.combined, state.details, std::array<int, 2>{ max_mm1 - state.first_match.second, max_mm2 - current2.second });

            if (state.details.index >= 0) {
                ++state.counts[state.details.index];
                return true;
            } else {
                return false;
            }
        };

        // Looping over all hits of the second for each hit of the first read.
        // This is done in a slightly convoluted way; we only search for
        // all hits of the second read _after_ we find the first hit of the
        // first read, so as to avoid a wasted search on the second read
        // if we never found a hit on the first read.
        while (inner_process(search_reverse1, constant1, max_mm1, against1.first, deets1, state.first_match)) {
            if (!deets2.finished) {
                // Alright, populating the second match buffer. We also
                // return immediately if any of them form a valid
                // combination with the first hit of the first read.
                while (inner_process(search_reverse2, constant2, max_mm2, against2.first, deets2, state.second_matches)) {
                    if (checker(state.second_matches.size() - 1)) {
                        return true;
                    }
                }
                if (state.second_matches.empty()) {
                    break;
                }
            } else {
                // And then this part does all the pairwise comparisons with
                // every hit in the first read.
                for (size_t i = 0; i < state.second_matches.size(); ++i) {
                    if (checker(i)) {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    std::pair<int, int> process_best(State& state, const std::pair<const char*, const char*>& against1, const std::pair<const char*, const char*>& against2) const {
        auto deets1 = constant1.initialize(against1.first, against1.second - against1.first);
        auto deets2 = constant2.initialize(against2.first, against2.second - against2.first);

        // Getting all hits on the second read, and then looping over that
        // vector for each hit of the first read. We have to do all pairwise
        // comparisons anyway to find the best hit.
        state.second_matches.clear();
        while (inner_process(search_reverse2, constant2, max_mm2, against2.first, deets2, state.second_matches)) {}

        int chosen = -1;
        int best_mismatches = max_mm1 + max_mm2 + 1;
        size_t num_second_matches = state.second_matches.size();

        if (!state.second_matches.empty()) {
            while (inner_process(search_reverse1, constant1, max_mm1, against1.first, deets1, state.first_match)) {
                for (size_t i = 0; i < num_second_matches; ++i) {
                    const auto& current2 = state.second_matches[i];

                    state.combined = state.first_match.first;
                    state.combined += current2.first; // separate line is deliberate.
                    varlib.search(state.combined, state.details, std::array<int, 2>{ max_mm1 - state.first_match.second, max_mm2 - current2.second });

                    if (state.details.index >= 0) {
                        int cur_mismatches = state.details.mismatches + state.first_match.second + current2.second;
                        if (cur_mismatches < best_mismatches) {
                            chosen = state.details.index;
                            best_mismatches = cur_mismatches;
                        } else if (cur_mismatches == best_mismatches && chosen != state.details.index) { // ambiguous.
                            chosen = -1;
                        }
                    }
                }
            }
        }

        return std::make_pair(chosen, best_mismatches);
    }

public:
    /**
     *@cond
     */
    bool process(State& state, const std::pair<const char*, const char*>& r1, const std::pair<const char*, const char*>& r2) const {
        bool found;

        if (use_first) {
            found = process_first(state, r1, r2);
            if (!found && randomized) {
                found = process_first(state, r2, r1);
            }

        } else {
            auto best = process_best(state, r1, r2);
            if (randomized) {
                auto best2 = process_best(state, r2, r1);
                if (best.first < 0 || best.second > best2.second) {
                    best = best2;
                } else if (best.second == best2.second && best.first != best2.first) {
                    best.first = -1; // ambiguous.
                }
            }

            found = best.first >= 0;
            if (found) {
                ++state.counts[best.first];
            }
        }

        ++state.total;
        return found;
    }
    /**
     *@endcond
     */

private:
    bool search_reverse1, search_reverse2;

    ScanTemplate<max_size> constant1, constant2;
    SegmentedBarcodeSearch<2> varlib;
    int max_mm1, max_mm2;

    bool randomized;
    bool use_first = true;

    std::vector<int> counts;
    int total = 0;

public:
    /**
     * @return Vector containing the frequency of each valid combination.
     * This has length equal to the number of valid dual barcode combinations (i.e., the length of `barcode_pool1` and `barcode_pool2` in the constructor).
     * Each entry contains the count for the corresponding dual barcode combination.
     */
    const std::vector<int>& get_counts() const {
        return counts;
    }

    /**
     * @return Total number of read pairs processed by the handler.
     */
    int get_total() const {
        return total;
    }
};

/**
 * @cond
 */
// Soft-deprecated back-compatible aliases.
template<size_t max_size>
using DualBarcodes = DualBarcodesPairedEnd<max_size>;
/**
 * @endcond
 */

}

#endif
