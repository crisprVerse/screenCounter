#ifndef BYTEME_OSTREAM_WRITER_HPP
#define BYTEME_OSTREAM_WRITER_HPP

#include <ostream>
#include <stdexcept>
#include "Writer.hpp"

/**
 * @file OstreamWriter.hpp
 *
 * @brief Write bytes to an arbitrary output stream.
 */

namespace byteme {

/**
 * @brief Read bytes from a `std::ostream`.
 *
 * This is just a wrapper around `std::ostream::write` for compatibility.
 */
class OstreamWriter : public Writer {
public:
    /**
     * @param output An output stream.
     * This is assumed to live until `finish()` is called.
     */
    OstreamWriter(std::ostream& output) : ptr(&output) {}

    void write(const unsigned char* buffer, size_t n) {
        ptr->write(reinterpret_cast<const char*>(buffer), n);
        if (!(ptr->good())) {
            throw std::runtime_error("failed to write to arbitrary output stream");
        }
    }

    void finish() {
        ptr->flush();
        if (ptr->fail() || ptr->bad()) {
            throw std::runtime_error("failed to flush to arbitrary output stream");
        }
    }

private:
    std::ostream* ptr;
};

}

#endif
