#ifndef KAORI_FASTQ_READER_HPP
#define KAORI_FASTQ_READER_HPP

#include "byteme/PerByte.hpp"
#include <cctype>
#include <vector>
#include <stdexcept>

/**
 * @file FastqReader.hpp
 *
 * @brief Defines the `FastqReader` class.
 */

namespace kaori {

/**
 * @brief Stream reads from a FASTQ file.
 *
 * Pretty much what it says on the tin.
 * Multi-line sequence and quality strings are supported.
 * The name of each read is only considered up to the first whitespace.
 */
class FastqReader {
public:
    /**
     * @param p Any `byteme::Reader` instance that defines a text stream.
     */
    FastqReader(byteme::Reader* p) : pb(p) {
        sequence.reserve(200);
        name.reserve(200);
        okay = pb.valid();
    }

    /**
     * Extract details for the next read in the file.
     *
     * @return Whether or not a record was successfully extracted.
     * If `true`, `get_sequence()` and `get_name()` may be used.
     * If `false`, this indicates that we reached the end of the file.
     */
    bool operator()() {
        // Quitting early if the buffer is already empty. 
        if (!okay) {
            return false;
        }

        int init_line = line_count;

        // Processing the name. This should be on a single line, hopefully.
        name.clear();
        char val = pb.get();
        if (val != '@') {
            throw std::runtime_error("read name should start with '@' (starting line " + std::to_string(init_line + 1) + ")");
        }

        val = advance_and_check();
        while (!std::isspace(val)) {
            name.push_back(val);
            val = advance_and_check();
        }

        while (val != '\n') {
            val = advance_and_check();
        }
        ++line_count;

        // Processing the sequence itself until we get to a '+'.
        sequence.clear();
        val = advance_and_check();
        while (val != '+') {
            if (val != '\n') {
                sequence.push_back(val);
            }
            val = advance_and_check();
        }
        ++line_count;

        // Line 3 should be a single line; starting with '+' is implicit from above.
        val = advance_and_check();
        while (val != '\n') {
            val = advance_and_check();
        } 
        ++line_count;

        // Processing the qualities. Extraction is allowed to fail if we're at
        // the end of the file. Note that we can't check for '@' as a
        // delimitor, as this can be a valid score, so instead we check at each
        // newline whether we've reached the specified length, and quit if so.
        size_t qual_length = 0, seq_length = sequence.size();
        okay = false;

        while (pb.advance()) {
            val = pb.get();
            if (val != '\n') {
                ++qual_length;
            } else if (qual_length >= seq_length) {
                okay = pb.advance(); // sneak past the newline.
                break;
            }
        }

        if (qual_length != seq_length) {
            throw std::runtime_error("non-equal lengths for quality and sequence strings (starting line " + std::to_string(init_line + 1) + ")");
        }

        ++line_count;

        return true;
    }

private:
    byteme::PerByte<> pb;

    char advance_and_check() {
        if (!pb.advance()) {
            throw std::runtime_error("premature end of the file at line " + std::to_string(line_count + 1));
        }
        return pb.get();
    }

private:
    std::vector<char> sequence;
    std::vector<char> name;
    bool okay;
    int line_count = 0;

public:
    /**
     * @return Vector containing the sequence for the current read.
     */
    const std::vector<char>& get_sequence() const {
        return sequence;
    }

    /**
     * @return Vector containing the name for the current read.
     * Note that the name is considered to end at the first whitespace on the line.
     */
    const std::vector<char>& get_name() const {
        return name;
    }
};

}

#endif
