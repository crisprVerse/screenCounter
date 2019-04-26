#include <algorithm>

void reverse_complement (char* seq, size_t len) {
    for (size_t i=0; i<len; ++i) {
        char& current=seq[i];
        switch(current) {
            case 'A': case 'a':
                current='T'; break;
            case 'C': case 'c':
                current='G'; break;
            case 'G': case 'g':
                current='C'; break;
            case 'T': case 't':
                current='A'; break;
        }
    }
    std::reverse(seq, seq+len);
    return;
}
