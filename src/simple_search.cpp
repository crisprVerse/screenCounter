#include "simple_search.h"
#include "hash_scanner.h"

combined_dictionary_part::const_iterator simple_search_internal(const combined_dictionary_part& ref, const char* seq, size_t seq_len, const size_t target_len) {
    combined_dictionary_part::const_iterator output = ref.end();

    if (seq_len >= target_len) {
        hash_scanner hasher(seq, target_len);
        size_t end = target_len;
        int best_nedits = target_len; // any large value will do.

        do {
            if (hasher.valid()) {
                auto it = ref.find(hasher.hash());
                if (it != ref.end()) {
                    int nedits = (it->second).first;
                    if (nedits < best_nedits) {
                        output = it;
                        best_nedits = nedits;
                        if (nedits == 0) {
                            break;
                        }
                    }
                }
            }

            if (end >= seq_len) {
                break;
            } 

            hasher.advance();
            ++end;
        } while (1);
    }

    return output;
}

void simple_search(const combined_dictionary& ref, const char* seq, size_t len, int& best_nedits, int& best_choice) {
    for (auto& ep : ref) {
        const auto& curdict = ep.second;
        auto selected = simple_search_internal(curdict, seq, len, ep.first); 
        if (selected != curdict.end() && (selected->second).first < best_nedits) {
            best_nedits = (selected->second).first;
            best_choice = (selected->second).second[0];
        } 
        if (best_nedits == 0) {
            break;
        }
    }
    return;
}

void simple_search(const combined_dictionary& ref, const char* seq, size_t len, int& best_nedits, std::vector<int>& best_combination) {
    for (const auto& ep : ref) {
        const auto& curdict = ep.second;
        auto selected = simple_search_internal(curdict, seq, len, ep.first); 
        if (selected != curdict.end() && (selected->second).first < best_nedits) {
            best_nedits = (selected->second).first;
            auto& combo = (selected->second).second;
            std::copy(combo.begin(), combo.end(), best_combination.begin());
        } 
        if (best_nedits == 0) {
            break;
        }
    }
    return;
}
