#ifndef KAORI_FASTQ_READER_HPP
#define KAORI_FASTQ_READER_HPP

#include "byteme/Reader.hpp"
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
    FastqReader(byteme::Reader* p) : ptr(p) {
        sequence.reserve(200);
        name.reserve(200);

        refresh();
        if (available) {
            if (buffer[0] != '@') {
                throw std::runtime_error("first line containing FASTQ name should start with '@'");
            }
            ++avail_pos;
        }
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
        if (available == 0) {
            return false;
        }
        int init_line = line_count;

        // Processing the name. This should be on a single line, hopefully.
        name.clear();
        while (1) {
            for (; avail_pos < available && !std::isspace(buffer[avail_pos]); ++avail_pos) {
                name.push_back(buffer[avail_pos]);
            }
            if (avail_pos == available) {
                refresh();
            } else {
                break;
            }
        }

        while (1) {
            bool finished = false;
            for (; avail_pos < available; ++avail_pos) {
                if (buffer[avail_pos] == '\n') {
                    finished = true;
                    break;
                }
            }
            if (!finished) {
                refresh();
            } else {
                ++line_count;
                ++avail_pos;
                break;
            }
        }

        // Processing the sequence itself.
        sequence.clear();
        while (1) {
            for (; avail_pos < available && buffer[avail_pos] != '\n'; ++avail_pos) {
                sequence.push_back(buffer[avail_pos]);
            }
            if (avail_pos == available) {
                refresh();
            } else {
                // Implicitly the current position is a newline.
                ++line_count;
                if (avail_pos + 1 < available) {
                    ++avail_pos;
                } else {
                    refresh();
                }
                if (buffer[avail_pos] == '+') {
                    break;
                }
            }
        }

        // Line 3 should be a single line. The check above implies that it
        // starts with '+', so no need to check it.
        while (1) {
            bool finished = false;
            for (; avail_pos < available; ++avail_pos) {
                if (buffer[avail_pos] == '\n') {
                    ++line_count;
                    finished = true;
                    break;
                }
            }
            if (!finished) {
                refresh();
            } else {
                ++line_count;
                ++avail_pos;
                break;
            }
        }

        // Processing the qualities. Extraction is allowed to fail if we're at the
        // end of the file.
        size_t qual_length = 0;

        while (1) {
            for (; avail_pos < available && buffer[avail_pos] != '\n'; ++avail_pos) {
                ++qual_length;
            }
            if (avail_pos == available) {
                refresh<false>();
                if (!available) {
                    break;
                }
            } else {
                // Implicitly the current position is a newline. 
                ++line_count;
                if (avail_pos + 1 < available) {
                    ++avail_pos;
                } else {
                    refresh<false>();
                    if (!available) {
                        break;
                    }
                }
                if (buffer[avail_pos] == '@') {
                    ++avail_pos; // skipping the marker for the next set of elements.
                    break;
                }
            }
        }

        if (qual_length != sequence.size()) {
            throw std::runtime_error("non-equal lengths for quality and sequence strings (starting line " + std::to_string(init_line) + ")");
        }

        return true;
    }

private:
    byteme::Reader* ptr;
    bool source_empty = false;

    const char * buffer;
    size_t available = 0;
    size_t avail_pos = 0;

    template<bool must_work = true>
    void refresh() {
        avail_pos = 0;

        if (source_empty) {
            if (must_work) {
                throw std::runtime_error("premature end of the file at line " + std::to_string(line_count + 1));
            } else {
                available = 0;
                return;
            }
        }

        source_empty = !(ptr->operator()());
        buffer = reinterpret_cast<const char*>(ptr->buffer());
        available = ptr->available();
    }

private:
    std::vector<char> sequence;
    std::vector<char> name;
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
