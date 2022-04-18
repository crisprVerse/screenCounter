#ifndef KAORI_SINGLE_BARCODE_SINGLE_END_HPP
#define KAORI_SINGLE_BARCODE_SINGLE_END_HPP

#include "../SimpleSingleMatch.hpp"
#include <vector>

namespace kaori {

template<size_t N>
class SingleBarcodeSingleEnd {
public:
    SingleBarcodeSingleEnd(const char* constant, size_t size, int strand, const std::vector<const char*>& variable, int mismatches = 0) : 
        matcher(constant, size, strand != 1, strand != 0, variable, mismatches), counts(variable.size()) {}
        
    SingleBarcodeSingleEnd& set_first(bool t = true) {
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
