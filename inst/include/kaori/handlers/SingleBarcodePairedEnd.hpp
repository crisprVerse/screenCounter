#ifndef KAORI_SINGLE_BARCODE_PAIRED_END_HPP
#define KAORI_SINGLE_BARCODE_PAIRED_END_HPP

#include "../SimpleSingleMatch.hpp"
#include <vector>

namespace kaori {

template<size_t N>
class SingleBarcodePairedEnd {
public:
    SingleBarcodePairedEnd(const char* constant, size_t size, int strand, const std::vector<const char*>& variable, int mismatches = 0) : 
        matcher(constant, size, strand != 1, strand != 0, variable, mismatches), counts(variable.size()) {}
        
    SingleBarcodePairedEnd& set_first(bool t = true) {
        use_first = t;
        return *this;
    }

public:
    struct State {
        State() {}

        State(typename SimpleSingleMatch<N>::SearchState s, size_t nvar) : search(std::move(s)), counts(nvar) {}

        typename SimpleSingleMatch<N>::SearchState search;
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

public:
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

private:
    SimpleSingleMatch<N> matcher;
    std::vector<int> counts;
    int total = 0;
    bool use_first = true;

public:
    const std::vector<int>& get_counts() const {
        return counts;        
    }

    int get_total() const {
        return total;
    }
};

}

#endif
