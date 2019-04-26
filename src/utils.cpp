#include <algorithm>

void reverse_complement (char* seq, size_t len) {
    for (size_t i=0; i<len; ++i) {
        char& current=seq[i];
        switch(current) {
            case 'A': case 'a':
                current='T'; return;
            case 'C': case 'c':
                current='G'; return;
            case 'G': case 'g':
                current='C'; return;
            case 'T': case 't':
                current='A'; return;
        }
    }
    std::reverse(seq, seq+len);
    return;
}
