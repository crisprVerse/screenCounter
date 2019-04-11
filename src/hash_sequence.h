#include <string>

std::u32string hash_sequence(const char*, const size_t);

void shift_sequence(std::u32string&, const size_t, char);

struct rolling_substitution {
public:
    rolling_substitution(const std::u32string&, size_t);
    bool advance();
    const std::u32string& get() const;
private:
    std::u32string hash;
    size_t len, word, pos;
    uint32_t current, original;
};

struct rolling_deletion {
public:
    rolling_deletion(const std::u32string&, size_t);
    bool advance();
    const std::u32string& get() const;
private:
    std::u32string hash;
    size_t len, word, pos;
    uint32_t current, discarded;
};
