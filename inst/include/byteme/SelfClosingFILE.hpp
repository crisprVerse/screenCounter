#ifndef BYTEME_SELF_CLOSING_FILE_HPP
#define BYTEME_SELF_CLOSING_FILE_HPP

#include <cstdio>
#include <string>
#include <stdexcept>

namespace byteme {

struct SelfClosingFILE {
    SelfClosingFILE(const char * path, const char* mode) : handle(std::fopen(path, mode)) {
        if (!handle) {
            throw std::runtime_error("failed to open file at '" + std::string(path) + "'");
        }
        return;
    }

    ~SelfClosingFILE() {
        if (handle) {
            // Don't check return value as we can't throw here anyway.
            std::fclose(handle);
        }
        return;
    }

    // Delete the remaining constructors.
    SelfClosingFILE(const SelfClosingFILE&) = delete;
    SelfClosingFILE(SelfClosingFILE&&) = delete;
    SelfClosingFILE& operator=(const SelfClosingFILE&) = delete;
    SelfClosingFILE& operator=(SelfClosingFILE&&) = delete;

    FILE* handle;
};

}

#endif
