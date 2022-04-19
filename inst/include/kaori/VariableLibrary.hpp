#ifndef KAORI_VARIABLE_LIBRARY_HPP
#define KAORI_VARIABLE_LIBRARY_HPP

#include "MismatchTrie.hpp"
#include <unordered_map>
#include <string>
#include <vector>

namespace kaori {

class VariableLibrary {
public:
    VariableLibrary() {}

    VariableLibrary(const std::vector<const char*>& options, size_t len, int mismatches = 0, bool reverse = false) : trie(len), max_mismatches(mismatches) {
        // The number of mismatches is an instance-specific constant to ensure
        // that the cache is valid across multiple calls to match(). The trie
        // search breaks early when it hits the cap, but the cap might be
        // different if the number of mismatches changes.  So if we break early
        // and report a miss, the miss will be cached and returned in cases
        // where there is a higher cap (and thus might actually be a hit).

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
                throw std::runtime_error("duplicate variable sequence '" + current + "'");
            }
            exact[current] = i;
            trie.add(current.c_str());
        }
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

public:
    void match(const std::string& x, SearchState& state) const {
        auto it = exact.find(x);
        if (it != exact.end()) {
            state.index = it->second;
            state.mismatches = 0;
            return;
        }

        // Seeing if it's any of the caches; otherwise searching the trie.
        std::pair<int, int> missed;
        auto cit = cache.find(x);

        if (cit == cache.end()) {
            auto lit = state.cache.find(x);
            if (lit != state.cache.end()) {
                missed = lit->second;
            } else {
                missed = trie.search(x.c_str(), max_mismatches);
                state.cache[x] = missed;
            }
        } else {
            missed = cit->second;
        }

        state.index = missed.first;
        state.mismatches = missed.second;
        return;
    }

private:
    std::unordered_map<std::string, int> exact;
    MismatchTrie trie;
    std::unordered_map<std::string, std::pair<int, int> > cache;
    int max_mismatches;
};

}

#endif
