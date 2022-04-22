#ifndef KAORI_COMBINATORIAL_BARCODES_PAIRED_END_HPP
#define KAORI_COMBINATORIAL_BARCODES_PAIRED_END_HPP

#include "../SimpleSingleMatch.hpp"
#include "../utils.hpp"

namespace kaori {

template<size_t N>
class CombinatorialBarcodesPairedEnd { 
public:
    CombinatorialBarcodesPairedEnd(
        const char* con1, size_t n1, bool rev1, const std::vector<const char*>& var1, int mm1, 
        const char* con2, size_t n2, bool rev2, const std::vector<const char*>& var2, int mm2,
        bool random = false,
        bool duplicates = false
    ) :
        matcher1(con1, n1, !rev1, rev1, var1, mm1, duplicates),
        matcher2(con2, n2, !rev2, rev2, var2, mm2, duplicates),
        randomized(random)
    {
        num_options[0] = var1.size();
        num_options[1] = var2.size();
    }

    CombinatorialBarcodesPairedEnd& set_first(bool t = true) {
        use_first = t;
        return *this;
    }

public:
    struct State {
        State() {}

        State(typename SimpleSingleMatch<N>::SearchState s1, typename SimpleSingleMatch<N>::SearchState s2) : search1(std::move(s1)), search2(std::move(s2)) {}

        std::vector<std::array<int, 2> >collected;
        int read1_only = 0;
        int read2_only = 0;
        int total = 0;

        /**
         * @cond
         */
        typename SimpleSingleMatch<N>::SearchState search1;
        typename SimpleSingleMatch<N>::SearchState search2;
        /**
         * @endcond
         */
    };

    State initialize() const {
        return State(matcher1.initialize(), matcher2.initialize());
    }

    void reduce(State& s) {
        matcher1.reduce(s.search1);
        matcher2.reduce(s.search2);
        combinations.insert(combinations.end(), s.collected.begin(), s.collected.end());
        total += s.total;
        read1_only += s.read1_only;
        read2_only += s.read2_only;
    }

    constexpr static bool use_names = false;

public:
    void process(State& state, const std::pair<const char*, const char*>& r1, const std::pair<const char*, const char*>& r2) const {
        if (use_first) {
            bool m1 = matcher1.search_first(r1.first, r1.second - r1.first, state.search1);
            bool m2 = matcher2.search_first(r2.first, r2.second - r2.first, state.search2);

            if (m1 && m2) {
                state.collected.emplace_back(std::array<int, 2>{ state.search1.index, state.search2.index });
            } else if (randomized) {
                bool n1 = matcher1.search_first(r2.first, r2.second - r2.first, state.search1);
                bool n2 = matcher2.search_first(r1.first, r1.second - r1.first, state.search2);
                if (n1 && n2) {
                    state.collected.emplace_back(std::array<int, 2>{ state.search1.index, state.search2.index });
                } else {
                    if (m1 || n1) {
                        ++state.read1_only;
                    } else if (m2 || n2) {
                        ++state.read2_only;
                    }
                }
            } else {
                if (m1) {
                    ++state.read1_only;
                } else if (m2) {
                    ++state.read2_only;
                }
            }

        } else {
            bool m1 = matcher1.search_best(r1.first, r1.second - r1.first, state.search1);
            bool m2 = matcher2.search_best(r2.first, r2.second - r2.first, state.search2);

            if (!randomized) {
                if (m1 && m2) {
                    state.collected.emplace_back(std::array<int, 2>{ state.search1.index, state.search2.index });
                } else if (m1) {
                    ++state.read1_only;
                } else if (m2) {
                    ++state.read2_only;
                }
            } else if (m1 && m2) {
                std::array<int, 2> candidate{ state.search1.index, state.search2.index };
                int mismatches = state.search1.mismatches + state.search2.mismatches;

                bool n1 = matcher1.search_best(r2.first, r2.second - r2.first, state.search1);
                bool n2 = matcher2.search_best(r1.first, r1.second - r1.first, state.search2);

                if (n1 && n2) {
                    int rmismatches = state.search1.mismatches + state.search2.mismatches;
                    if (mismatches > rmismatches) {
                        state.collected.emplace_back(std::array<int, 2>{ state.search1.index, state.search2.index });
                    } else if (mismatches < rmismatches) {
                        state.collected.emplace_back(candidate);
                    } else if (candidate[0] == state.search1.index && candidate[1] == state.search2.index) {
                        // If the mismatches are the same, it may not be ambiguous
                        // if the indices would be the same anyway.
                        state.collected.emplace_back(candidate);
                    }
                } else {
                    state.collected.emplace_back(candidate);
                }
            } else {
                bool n1 = matcher1.search_best(r2.first, r2.second - r2.first, state.search1);
                bool n2 = matcher2.search_best(r1.first, r1.second - r1.first, state.search2);

                if (n1 && n2) {
                    state.collected.emplace_back(std::array<int, 2>{ state.search1.index, state.search2.index });
                } else if (m1 || n1) {
                    ++state.read1_only;
                } else if (m2 || n2) {
                    ++state.read2_only;
                }
            }
        }

        ++state.total;
    }

public:
    void sort() {
        sort_combinations(combinations, num_options);
    }

    const std::vector<std::array<int, 2> >& get_combinations() const {
        return combinations;
    }

    int get_total() const {
        return total;
    }

    int get_read1_only() const {
        return read1_only;
    }

    int get_read2_only() const {
        return read2_only;
    }
private:
    SimpleSingleMatch<N> matcher1, matcher2;
    std::array<size_t, 2> num_options;

    bool randomized;
    bool use_first = true;

    std::vector<std::array<int, 2> > combinations;
    int total = 0;
    int read1_only = 0;
    int read2_only = 0;
};

}

#endif
