#include <string>

std::u32string hash_sequence(const char*, const size_t);

void shift_sequence(std::u32string&, const size_t, char);

bool is_valid (char);

int is_valid (const char*, size_t);

class hash_scanner {
public:
    hash_scanner(const char*, size_t);
    void advance ();
    const std::u32string& hash() const;
    bool valid() const;
private:
    const char* ptr;
    size_t len;
    std::u32string hashed;
    size_t nvalid;
};

