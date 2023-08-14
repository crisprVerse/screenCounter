#ifndef KAORI_BARCODE_POOL_HPP
#define KAORI_BARCODE_POOL_HPP

#include <vector>
#include <string>

/**
 * @file BarcodePool.hpp
 *
 * @brief Defines the `BarcodePool` class.
 */

namespace kaori {

/**
 * @brief Pool of barcode sequences for a variable region.
 *
 * The `BarcodePool` class defines the pool of possible barcode sequences for a given variable region in the template sequence.
 * All sequences in this set are assumed to have the same length.
 */
struct BarcodePool {
    /**
     * Default constructor.
     */
    BarcodePool() {}
    
    /**
     * @param barcode_pool Vector of pointers to sequences of length `l`, containing the pool of possible barcode sequences.
     * @param barcode_length Length of each sequence.
     */
    BarcodePool(std::vector<const char*> barcode_pool, size_t barcode_length) : pool(std::move(barcode_pool)), length(barcode_length) {}

    /**
     * @param barcode_pool Vector of sequences of the same length, containing the pool of possible barcode sequences.
     *
     * It is assumed that the lifetime of `barcode_pool` (and its strings) exceeds that of the constructed `BarcodePool`. 
     */
    BarcodePool(const std::vector<std::string>& barcode_pool) {
        if (barcode_pool.size()) {
            length = barcode_pool.front().size();
            pool.reserve(barcode_pool.size());
            for (const auto& x : barcode_pool) {
                if (x.size() != length) {
                    throw std::runtime_error("sequences for a given variable region should be of a constant length");
                }
                pool.push_back(x.c_str());
            }
        }
    }

    /**
     * Vector containing pointers to sequences of length `len`.
     */
    std::vector<const char*> pool;

    /**
     * Length of each sequence in `pool`.
     */
    size_t length = 0;

    /**
     * @return Number of sequences in the pool.
     */
    size_t size() const {
        return pool.size();
    }

    /**
     * @param i Index of the barcode  sequence of interest.
     * @return Pointer to the `i`-th sequence in the pool.
     */
    const char* operator[](size_t i) const {
        return pool[i];
    }
};

}

#endif
