#ifndef HASH_SCANNER_H
#define HASH_SCANNER_H

#include "hash_sequence.h"

class hash_scanner {
public:
    hash_scanner() = default;
    hash_scanner(const char*, size_t);
    void advance ();
    const seqhash& hash() const;
    bool valid() const;
private:
    const char* ptr = NULL;
    size_t len = 0;
    seqhash hashed;
    size_t nvalid = 0;
};

#endif
