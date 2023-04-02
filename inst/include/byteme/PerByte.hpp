#ifndef BYTEME_PERBYTE_HPP
#define BYTEME_PERBYTE_HPP

/**
 * @file PerByte.hpp
 * @brief Perform byte-by-byte extraction.
 */

#include <thread>
#include <vector>
#include <algorithm>
#include <exception>

#include "Reader.hpp"

namespace byteme {

/**
 * @brief Perform byte-by-byte extraction from a `Reader` source.
 *
 * @tparam T Character type to return, usually `char` for text or `unsigned char` for binary.
 *
 * This wraps a `Reader` so that developers can avoid the boilerplate of managing blocks of bytes,
 * when all they want is to iterate over the bytes of the input.
 */
template<typename T = char>
struct PerByte {
private:
    const T* ptr = nullptr;
    size_t available = 0;
    size_t current = 0;
    bool remaining = false;
    size_t overall = 0;

    Reader* reader = nullptr;

    void refill() {
        remaining = (*reader)();
        ptr = reinterpret_cast<const T*>(reader->buffer());
        available = reader->available();
        current = 0;
    }
public:
    /**
     * @param r An existing reader object that has not been read from.
     */
    PerByte(Reader& r) : reader(&r) {
        refill();
    }

    /**
     * @return Whether this instance still has bytes to be read.
     */
    bool valid() const {
        return current < available;
    }

    /**
     * Advance to the next byte, possibly reading a new chunk from the supplied `Reader`.
     * This should only be called if `valid()` is `true`.
     */
    void advance() {
        ++current;
        if (current == available) {
            overall += available;
            if (remaining) {
                refill();
            }
        }
    }

    /**
     * @return The current byte.
     *
     * This should only be called if `valid()` is `true`.
     */
    T get() const {
        return ptr[current];
    }

    /**
     * @return The position of the current byte since the start of the input.
     */
    size_t position() const {
        return overall + current;
    }
};

/**
 * @brief Perform parallelized byte-by-byte extraction from a `Reader` source.
 *
 * @tparam T Character type to return, usually `char` for text or `unsigned char` for binary.
 *
 * This is much like `PerByte` except that the `Reader`'s loading operation is called in a separate thread,
 * thus allowing the caller to parse the bytes of the current chunk in parallel.
 */
template<typename T = char>
struct PerByteParallel {
private:
    size_t current = 0;
    bool remaining = false;
    size_t available = 0;
    size_t overall = 0;

    Reader* reader = nullptr;
    bool use_meanwhile = false;
    std::thread meanwhile;
    std::exception_ptr thread_err = nullptr;
    std::vector<T> buffer;

    void refill() {
        auto ptr = reinterpret_cast<const T*>(reader->buffer());
        available = reader->available();
        buffer.resize(available);
        std::copy(ptr, ptr + available, buffer.begin());

        current = 0;
        use_meanwhile = remaining;
        if (remaining) {
            meanwhile = std::thread([&]() -> void {
                try {
                    remaining = (*reader)(); 
                } catch (...) {
                    thread_err = std::current_exception();
                }
            });
        }
    }

public:
    /**
     * @param r An existing reader object that has not been read from.
     */
    PerByteParallel(Reader& r) : reader(&r) {
        remaining = (*reader)();
        refill();
    }

    /**
     * @cond
     */
    ~PerByteParallel() {
        // Avoid a dangling thread if we need to destroy this prematurely,
        // e.g., because the caller encountered an exception.
        if (use_meanwhile) {
            meanwhile.join();
        }
    }
    /**
     * @endcond
     */

    /**
     * @return Whether this instance still has bytes to be read.
     */
    bool valid() const {
        return current < available;
    }

    /**
     * Advance to the next byte, possibly reading a new chunk from the supplied `Reader`.
     * This should only be called if `valid()` is `true`.
     */
    void advance() {
        ++current;
        if (current == available) {
            overall += available;
            if (use_meanwhile) {
                meanwhile.join();
                if (thread_err) {
                    std::rethrow_exception(thread_err);
                }
                refill();
            }
        }
    }

    /**
     * @return The current byte.
     *
     * This should only be called if `valid()` is `true`.
     */
    T get() const {
        return buffer[current];
    }

    /**
     * @return The position of the current byte since the start of the input.
     */
    size_t position() const {
        return overall + current;
    }
};

}

#endif
