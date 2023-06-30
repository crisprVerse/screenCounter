#ifndef KAORI_DUAL_BARCODES_HPP
#define KAORI_DUAL_BARCODES_HPP

#include "../ScanTemplate.hpp"
#include "../BarcodeSearch.hpp"
#include "../utils.hpp"

/**
 * @file DualBarcodes.hpp
 *
 * @brief Process dual barcodes.
 */

namespace kaori {

/**
 * @brief Handler for dual barcodes.
 *
 * In this design, each read contains a target sequence created from a template with a single variable region.
 * For one read, the barcode is drawn from one pool of options, while the other read contains a barcode from another pool.
 * However, unlike `CombinatorialBarcodesPairedEnd`, the combinations are not random but are specifically assembled, typically corresponding to specific pairs of genes.
 * This handler will capture the frequencies of each barcode combination. 
 *
 * @tparam max_size Maximum length of the template sequences on both reads.
 */
template<size_t max_size>
class DualBarcodes { 
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
    DualBarcodes(
        const char* template_seq1, size_t template_length1, bool reverse1, const BarcodePool& barcode_pool1, int max_mismatches1, 
        const char* template_seq2, size_t template_length2, bool reverse2, const BarcodePool& barcode_pool2, int max_mismatches2,
        bool random = false
    ) :
        search_reverse1(reverse1),
        search_reverse2(reverse2),
        constant1(template_seq1, template_length1, !search_reverse1, search_reverse1),
        constant2(template_seq2, template_length2, !search_reverse2, search_reverse2),
        max_mm1(max_mismatches1),
        max_mm2(max_mismatches2),
        randomized(random)
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
                    current += reverse_complement(ptr1[len1 - j - 1]);
                }
            } else {
                current.insert(current.end(), ptr1, ptr1 + len1);
            }

            auto ptr2 = barcode_pool2[i];
            if (search_reverse2) {
                for (size_t j = 0; j < len2; ++j) {
                    current += reverse_complement(ptr2[len2 - j - 1]);
                }
            } else {
                current.insert(current.end(), ptr2, ptr2 + len2);
            }

            combined.push_back(std::move(current));
        }

        // Constructing the combined varlib.
        BarcodePool combined_set(combined);
        varlib = SegmentedBarcodeSearch(
            combined_set,
            std::array<int, 2>{ static_cast<int>(len1), static_cast<int>(len2) }, 
            std::array<int, 2>{ max_mm1, max_mm2 }
        );
    }

    /**
     * @param t Whether to search only for the first match to valid target sequence(s) across both reads.
     * If `false`, the handler will search for the best match (i.e., fewest mismatches) instead.
     *
     * @return A reference to this `DualBarcodes` instance.
     */
    DualBarcodes& set_first(bool t = true) {
        use_first = t;
        return *this;
    }

public:
    /**
     *@cond
     */
    struct State {
        State(size_t n = 0) : counts(n) {}
        std::vector<int> counts;
        int total = 0;

        std::vector<std::pair<std::string, int> > buffer2;

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
    static void emit_output(std::pair<std::string, int>& output, const char* start, const char* end, int mm) {
        output.first = std::string(start, end);
        output.second = mm;
        return;
    }

    static void emit_output(std::vector<std::pair<std::string, int> >& output, const char* start, const char* end, int mm) {
        output.emplace_back(std::string(start, end), mm);
        return;
    }

    template<class Store>
    static bool inner_process(
        bool reverse, 
        const ScanTemplate<max_size>& constant, 
        int max_mm,
        const char* against,
        typename ScanTemplate<max_size>::State& deets,
        Store& output)
    {
        while (!deets.finished) {
            constant.next(deets);
            if (reverse) {
                if (deets.reverse_mismatches <= max_mm) {
                    const auto& reg = constant.template variable_regions<true>()[0];
                    auto start = against + deets.position;
                    emit_output(output, start + reg.first, start + reg.second, deets.reverse_mismatches);
                    return true;
                }
            } else {
                if (deets.forward_mismatches <= max_mm) {
                    const auto& reg = constant.variable_regions()[0];
                    auto start = against + deets.position;
                    emit_output(output, start + reg.first, start + reg.second, deets.forward_mismatches);
                    return true;
                }
            }
        }
        return false;
    }

    bool process_first(State& state, const std::pair<const char*, const char*>& against1, const std::pair<const char*, const char*>& against2) const {
        auto deets1 = constant1.initialize(against1.first, against1.second - against1.first);
        std::pair<std::string, int> match1;

        auto deets2 = constant2.initialize(against2.first, against2.second - against2.first);
        state.buffer2.clear();

        auto checker = [&](size_t idx2) -> bool {
            const auto& current2 = state.buffer2[idx2];
            auto combined = match1.first + current2.first;
            varlib.search(combined, state.details, std::array<int, 2>{ max_mm1 - match1.second, max_mm2 - current2.second });

            if (state.details.index != -1) {
                ++state.counts[state.details.index];
                return true;
            } else {
                return false;
            }
        };

        // Looping over all hits of the second for each hit of the first.
        while (inner_process(search_reverse1, constant1, max_mm1, against1.first, deets1, match1)) {
            if (deets2.finished) {
                for (size_t i = 0; i < state.buffer2.size(); ++i) {
                    if (checker(i)) {
                        return true;
                    }
                }
            } else {
                while (inner_process(search_reverse2, constant2, max_mm2, against2.first, deets2, state.buffer2)) {
                    if (checker(state.buffer2.size() - 1)) {
                        return true;
                    }
                }
                if (state.buffer2.empty()) {
                    break;
                }
            }
        }

        return false;
    }

    std::pair<int, int> process_best(State& state, const std::pair<const char*, const char*>& against1, const std::pair<const char*, const char*>& against2) const {
        auto deets1 = constant1.initialize(against1.first, against1.second - against1.first);
        std::pair<std::string, int> match1;

        auto deets2 = constant2.initialize(against2.first, against2.second - against2.first);
        state.buffer2.clear();

        int chosen = -1;
        int best_mismatches = max_mm1 + max_mm2 + 1;

        auto checker = [&](size_t idx2) -> void {
            const auto& current2 = state.buffer2[idx2];
            auto combined = match1.first + current2.first;
            varlib.search(combined, state.details, std::array<int, 2>{ max_mm1 - match1.second, max_mm2 - current2.second });

            int cur_mismatches = state.details.mismatches;
            if (cur_mismatches < best_mismatches) {
                chosen = state.details.index;
                best_mismatches = cur_mismatches;
            } else if (cur_mismatches == best_mismatches && chosen != state.details.index) { // ambiguous.
                chosen = -1;
            }
        };

        while (inner_process(search_reverse1, constant1, max_mm1, against1.first, deets1, match1)) {
            if (deets2.finished) {
                for (size_t i = 0; i < state.buffer2.size(); ++i) {
                    checker(i);
                }
            } else {
                while (inner_process(search_reverse2, constant2, max_mm2, against2.first, deets2, state.buffer2)) {
                    checker(state.buffer2.size() - 1);
                }
                if (state.buffer2.empty()) {
                    break;
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

}

#endif
