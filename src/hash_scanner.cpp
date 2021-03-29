#include "hash_scanner.h"
#include "utils.h"

/* Class that combines shift_sequence and is_valid to iterate across a
 * string and produce the hash at each position (or die trying).
 */

hash_scanner::hash_scanner(const char* p, size_t n) : ptr(p), len(n), 
    hashed(hash_sequence(p, n)), nvalid(is_valid(p, n)) {}

void hash_scanner::advance() {
    nvalid-=is_valid(*ptr);
    const char next=*(ptr+len);
    ++ptr;

    shift_sequence(hashed, len, next);
    nvalid+=is_valid(next);
    return;
}

const seqhash& hash_scanner::hash() const { return hashed; }

bool hash_scanner::valid() const { return nvalid==len; }
