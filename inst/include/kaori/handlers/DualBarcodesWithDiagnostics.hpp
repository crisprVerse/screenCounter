#ifndef KAORI_DUAL_BARCODES_WITH_DIAGNOSTICS_HPP
#define KAORI_DUAL_BARCODES_WITH_DIAGNOSTICS_HPP

#include "DualBarcodes.hpp"
#include "CombinatorialBarcodesPairedEnd.hpp"
#include "../utils.hpp"

namespace kaori {

template<size_t N>
class DualBarcodesWithDiagnostics { 
public:
    DualBarcodesWithDiagnostics(
        const char* con1, size_t n1, bool rev1, const std::vector<const char*>& var1, int mm1, 
        const char* con2, size_t n2, bool rev2, const std::vector<const char*>& var2, int mm2,
        bool random = false
    ) :
        dual_handler(con1, n1, rev1, var1, mm1, con2, n2, rev2, var2, mm2, random),
        combo_handler(con1, n1, rev1, var1, mm1, con2, n2, rev2, var2, mm2, random, true) // we allow duplicates in the trie.
    {}

    DualBarcodesWithDiagnostics& set_first(bool t = true) {
        dual_handler.set_first(t);
        combo_handler.set_first(t);
        return *this;
    }

private:
    DualBarcodes<N> dual_handler;
    CombinatorialBarcodesPairedEnd<N> combo_handler;

public:
    struct State {
        State() {}
        State(typename DualBarcodes<N>::State ds, typename CombinatorialBarcodesPairedEnd<N>::State cs) : dual_state(std::move(ds)), combo_state(std::move(cs)) {}

        /**
         * @cond
         */
        typename DualBarcodes<N>::State dual_state;
        typename CombinatorialBarcodesPairedEnd<N>::State combo_state;
        /**
         * @endcond
         */
    };

    State initialize() const {
        return State(dual_handler.initialize(), combo_handler.initialize());
    }

    void reduce(State& s) {
        dual_handler.reduce(s.dual_state);
        combo_handler.reduce(s.combo_state);
    }

    constexpr static bool use_names = false;

public:
    void process(State& state, const std::pair<const char*, const char*>& r1, const std::pair<const char*, const char*>& r2) const {
        // Only searching for combinations if we couldn't find a proper dual barcode match.
        if (!dual_handler.process(state.dual_state, r1, r2)) {
            combo_handler.process(state.combo_state, r1, r2);
        }
    }

public:
    void sort() {
        combo_handler.sort();
    }

    const std::vector<int>& get_counts() const {
        return dual_handler.get_counts();
    }

    const std::vector<std::array<int, 2> >& get_combinations() const {
        return combo_handler.get_combinations();
    }

    int get_total() const {
        return dual_handler.get_total();
    }

    int get_read1_only() const {
        return combo_handler.get_read1_only();
    }

    int get_read2_only() const {
        return combo_handler.get_read2_only();
    }
};

}

#endif
