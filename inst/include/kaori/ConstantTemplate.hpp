#ifndef KAORI_CONSTANT_TEMPLATE_HPP
#define KAORI_CONSTANT_TEMPLATE_HPP

#include <bitset>
#include <deque>
#include "utils.hpp"

namespace kaori {

template<size_t N>
class ConstantTemplate { 
public:
    ConstantTemplate() {}

    ConstantTemplate(const char* s, size_t n, bool f, bool r) : length(n), forward(f), reverse(r) {
        if (n * 4 > N) {
            throw std::runtime_error("maximum constant size should be " + std::to_string(N/4) + " nt");
        }

        if (forward) {
            for (size_t i = 0; i < n; ++i) {
                char b = s[i];
                if (b != '-') {
                    add_base(forward_ref, b);
                    add_mask(forward_mask, i);
                } else {
                    shift(forward_ref);
                    shift(forward_mask);
                    add_variable_base(forward_variables, i);
                }
            }
        } else {
            // Forward variable regions are always defined.
            for (size_t i = 0; i < n; ++i) {
                char b = s[i];
                if (b == '-') {
                    add_variable_base(forward_variables, i);
                }
            }
        }

        if (reverse) {
            for (size_t i = 0; i < n; ++i) {
                char b = s[n - i - 1];
                if (b != '-') {
                    add_base(reverse_ref, reverse_complement(b));
                    add_mask(reverse_mask, i);
                } else {
                    shift(reverse_ref);
                    shift(reverse_mask);
                    add_variable_base(reverse_variables, i);
                }
            }
        }
    }

public:
    struct MatchDetails {
        size_t position = static_cast<size_t>(-1); // overflow should be sane.
        int forward_mismatches = -1;
        int reverse_mismatches = -1;
        bool finished = false;

        /**
         * @cond
         */
        std::bitset<N> state, ambiguous;
        const char * seq;
        size_t len;
        std::deque<size_t> bad;
        /**
         * @endcond
         */
    };

    MatchDetails initialize(const char* seq, size_t len) const {
        MatchDetails out;
        out.seq = seq;
        out.len = len;

        if (length <= len) {
            for (size_t i = 0; i < length - 1; ++i) {
                char base = seq[i];

                if (is_good(base)) {
                    add_base(out.state, base);
                    if (!out.bad.empty()) {
                        shift(out.ambiguous);
                    }
                } else {
                    shift(out.state);
                    out.state |= other_<N>;
                    shift(out.ambiguous);
                    out.ambiguous |= other_<N>;
                    out.bad.push_back(i);
                }
            }
        } else {
            out.finished = true;
        }

        return out;
    }

    void next(MatchDetails& match) const {
        if (!match.bad.empty() && match.bad.front() == match.position) {
            match.bad.pop_front();
            if (match.bad.empty()) {
                // This should effectively clear the ambiguous bitset, allowing
                // us to skip its shifting if there are no more ambiguous
                // bases. We do it here because we won't get an opportunity to
                // do it later; as 'bad' is empty, the shift below is skipped.
                shift(match.ambiguous); 
            }
        }

        size_t right = match.position + length;
        char base = match.seq[right];
        if (is_good(base)) {
            add_base(match.state, base); // no need to trim off the end, the mask will handle that.
            if (!match.bad.empty()) {
                shift(match.ambiguous);
            }
        } else {
            shift(match.state);
            match.state |= other_<N>;
            shift(match.ambiguous);
            match.ambiguous |= other_<N>;
            match.bad.push_back(right);
        }

        ++match.position;
        full_match(match);
        if (right + 1 == match.len) {
            match.finished = true;
        }

        return;
    }

private:
    std::bitset<N> forward_ref, forward_mask;
    std::bitset<N> reverse_ref, reverse_mask;
    size_t length;
    int mismatches;
    bool forward, reverse;

    static void add_mask(std::bitset<N>& current, size_t pos) {
        shift(current);
        current |= other_<N>;
        for (int i = 0; i < 4; ++i) {
            current[i] = 1;
        }
    }

    static int strand_match(const MatchDetails& match, const std::bitset<N>& ref, const std::bitset<N>& mask) {
        // pop count here is equal to the number of non-ambiguous mismatches *
        // 2 + number of ambiguous mismatches * 3. This is because
        // non-ambiguous bases are encoded by 1 set bit per 4 bases (so 2 are
        // left after a XOR'd mismatch), while ambiguous mismatches are encoded
        // by all set bits per 4 bases (which means that 3 are left after XOR).
        int pcount = ((match.state & mask) ^ ref).count(); 

        // Counting the number of ambiguous bases after masking. Each ambiguous
        // base is represented by 4 set bits, so we divide by 4 to get the number
        // of bases; then we multiply by three to remove their contribution. The
        // difference is then divided by two to get the number of non-ambig mm's.
        if (!match.bad.empty()) {
            int acount = (match.ambiguous & mask).count();
            acount /= 4;
            return acount + (pcount - acount * 3) / 2;
        } else {
            return pcount / 2;
        }
    }

    void full_match(MatchDetails& match) const {
        if (forward) {
            match.forward_mismatches = strand_match(match, forward_ref, forward_mask);
        }
        if (reverse) {
            match.reverse_mismatches = strand_match(match, reverse_ref, reverse_mask);
        }
    }

public:
    std::vector<std::pair<int, int> > forward_variables, reverse_variables;

    const std::vector<std::pair<int, int> >& variable_regions(bool reverse = false) const {
        return (reverse ? reverse_variables : forward_variables);
    } 

private:
    static void add_variable_base(std::vector<std::pair<int, int> >& variables, int i) {
        if (!variables.empty()) {
            auto& last = variables.back().second;
            if (last == i) {
                ++last;
                return;
            }
        }
        variables.emplace_back(i, i + 1);
        return;
    }
};

}

#endif
