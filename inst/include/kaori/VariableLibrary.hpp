#ifndef KAORI_VARIABLE_LIBRARY_HPP
#define KAORI_VARIABLE_LIBRARY_HPP

#include "MismatchTrie.hpp"
#include "utils.hpp"
#include <unordered_map>
#include <string>
#include <vector>
#include <array>

namespace kaori {

/** 
 * @cond
 */
template<class Trie>
void fill_library(
    const std::vector<const char*>& options, 
    std::unordered_map<std::string, int>& exact,
    Trie& trie,
    bool reverse,
    bool duplicates
) {
    size_t len = trie.get_length();

    for (size_t i = 0; i < options.size(); ++i) {
        auto ptr = options[i];

        std::string current;
        if (!reverse) {
            current = std::string(ptr, ptr + len);
        } else {
            for (int j = 0; j < len; ++j) {
                current += reverse_complement(ptr[len - j - 1]);
            }
        }

        auto it = exact.find(current);
        if (exact.find(current) != exact.end()) {
            if (!duplicates) {
                throw std::runtime_error("duplicate variable sequence '" + current + "'");
            }
        } else {
            exact[current] = i;
        }

        // Note that this must be called, even if the sequence is duplicated;
        // otherwise the trie's internal counter will not be properly incremented.
        trie.add(current.c_str(), duplicates);
    }
    return;
}

template<class Indexer, class Updater, class Cache, class Trie, class Result, class Mismatch>
void matcher_in_the_rye(const std::string& x, const Cache& cache, const Trie& trie, Result& res, const Mismatch& mismatches, const Mismatch& max_mismatches) {
    // Seeing if it's any of the caches; otherwise searching the trie.
    auto cit = cache.find(x);
    if (cit == cache.end()) {
        auto lit = res.cache.find(x);
        if (lit != res.cache.end()) {
            Updater::update(res, lit->second);
        } else {
            auto missed = trie.search(x.c_str(), mismatches);

            // The trie search breaks early when it hits the mismatch cap,
            // but the cap might be different across calls. If we break
            // early and report a miss, the miss will be cached and
            // returned in cases where there is a higher cap (and thus
            // might actually be a hit). As such, we should only store a
            // miss in the cache when the requested number of mismatches is
            // equal to the maximum value specified in the constructor.
            if (Indexer::index(missed) >= 0 || mismatches == max_mismatches) {
                res.cache[x] = missed;
            }

            Updater::update(res, missed);
        }
    } else {
        Updater::update(res, cit->second);
    }
    return;
}
/** 
 * @endcond
 */

class SimpleVariableLibrary {
public:
    SimpleVariableLibrary() {}

    SimpleVariableLibrary(const std::vector<const char*>& options, size_t len, int m = 0, bool reverse = false, bool duplicates = false) : trie(len), max_mismatches(m) {
        fill_library(options, exact, trie, reverse, duplicates);
        return;
    }

public:
    struct SearchState {
        int index = 0;
        int mismatches = 0;
        
        /**
         * @cond
         */
        std::unordered_map<std::string, std::pair<int, int> > cache;
        /**
         * @endcond
         */
    };

    SearchState initialize() const {
        return SearchState();
    }

    void reduce(SearchState& state) {
        cache.merge(state.cache);
        state.cache.clear();
    }

private:
    struct Index {
        static int index(const std::pair<int, int>& val) {
            return val.first;
        }
    };

    struct Updator {
        static void update(SearchState& state, const std::pair<int, int>& val) {
            state.index = val.first;
            state.mismatches = val.second;
            return;
        }
    };

public:
    void match(const std::string& x, SearchState& state) const {
        match(x, state, max_mismatches);
        return;
    }

    void match(const std::string& x, SearchState& state, int mismatches) const {
        auto it = exact.find(x);
        if (it != exact.end()) {
            state.index = it->second;
            state.mismatches = 0;
        } else {
            matcher_in_the_rye<Index, Updator>(x, cache, trie, state, mismatches, max_mismatches);
        }
    }

private:
    std::unordered_map<std::string, int> exact;
    SimpleMismatchTrie trie;
    std::unordered_map<std::string, std::pair<int, int> > cache;
    int max_mismatches;
};

template<size_t num_segments>
class SegmentedVariableLibrary {
public:
    SegmentedVariableLibrary() {}

    SegmentedVariableLibrary(const std::vector<const char*>& options, std::array<int, num_segments> segments, std::array<int, num_segments> mismatches, bool reverse = false, bool duplicates = false) : 
        trie(segments), 
        max_mismatches(mismatches) 
    {
        fill_library(options, exact, trie, reverse, duplicates);
        return;
    }

public:
    struct SearchState {
        SearchState() : per_segment() {}
        int index = 0;
        int mismatches = 0;
        std::array<int, num_segments> per_segment;
        
        /**
         * @cond
         */
        std::unordered_map<std::string, typename SegmentedMismatchTrie<num_segments>::SearchResult> cache;
        /**
         * @endcond
         */
    };

    SearchState initialize() const {
        return SearchState();
    }

    void reduce(SearchState& state) {
        cache.merge(state.cache);
        state.cache.clear();
    }

private:
    typedef typename SegmentedMismatchTrie<num_segments>::SearchResult SegmentedResult;

    struct Index {
        static int index(const SegmentedResult& val) {
            return val.index;
        }
    };

    struct Updator {
        static void update(SearchState& state, const SegmentedResult& val) {
            state.index = val.index;
            state.mismatches = val.total;
            state.per_segment = val.per_segment;
            return;
        }
    };

public:
    void match(const std::string& x, SearchState& state) const {
        match(x, state, max_mismatches);
        return;
    }

    void match(const std::string& x, SearchState& state, std::array<int, num_segments> mismatches) const {
        auto it = exact.find(x);
        if (it != exact.end()) {
            state.index = it->second;
            state.mismatches = 0;
            std::fill_n(state.per_segment.begin(), num_segments, 0);
        } else {
            matcher_in_the_rye<Index, Updator>(x, cache, trie, state, mismatches, max_mismatches);
        }
    }

private:
    std::unordered_map<std::string, int> exact;
    SegmentedMismatchTrie<num_segments> trie;
    std::unordered_map<std::string, SegmentedResult> cache;
    std::array<int, num_segments> max_mismatches;
};

}

#endif
