#ifndef SEARCH_SEQUENCE_H
#define SEARCH_SEQUENCE_H

#include <string>
#include "build_dictionary.h"
#include "hash_sequence.h"

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

/* Sets up search information. */

template<size_t nvariable>
struct combination {
    int indices[nvariable];
};

template<size_t nvariable>
struct search_info {
    search_info(Rcpp::List, Rcpp::StringVector, Rcpp::LogicalVector, Rcpp::LogicalVector);
    std::vector<sequence_dictionary> perfect, deleted;
    std::vector<std::u32string> constant_hash;
    std::vector<int> constant_lengths, variable_lengths;
    std::vector<size_t> constant_starts, variable_starts;
    size_t total_len;
    bool allow_del;
};

template<size_t nvariable>
search_info<nvariable>::search_info(Rcpp::List Guides, Rcpp::StringVector Constants,
    Rcpp::LogicalVector allow_subs, Rcpp::LogicalVector allow_dels) 
{
    if (allow_subs.size()!=1) {
        throw std::runtime_error("should be a logical scalar");
    }
    const bool allow_sub=allow_subs[0];
    allow_del=allow_dels[0];

    // Setting up the guides.
    if (nvariable!=Guides.size()) {
        std::stringstream err;
        err << "expecting " << nvariable << " variable regions";
        throw std::runtime_error(err.str());
    }
    perfect.reserve(nvariable);
    deleted.reserve(nvariable);
    variable_lengths.reserve(nvariable);

    for (size_t g=0; g<nvariable; ++g) {
        Rcpp::StringVector guides(Guides[g]);
        perfect.push_back(build_dictionary(guides, allow_sub));
        variable_lengths.push_back(perfect.back().len);
        if (allow_del) {
            deleted.push_back(build_deleted_dictionary(guides));
        }
    }

    // Setting up the constant regions.
    constexpr size_t nconstant=nvariable+1;
    if (Constants.size()!=nconstant) {
        throw std::runtime_error("number of constant regions should be 1 more than variable regions");
    }

    constant_hash.reserve(nconstant);
    constant_lengths.reserve(nconstant);

    for (size_t s=0; s<nconstant; ++s) {
        Rcpp::String con(Constants[s]);
        const size_t len=Rf_length(con.get_sexp());
        const char* ptr=con.get_cstring();
        constant_hash.push_back(hash_sequence(ptr, len)); 
        constant_lengths.push_back(len);
    }

    // Setting up search parameters.
    constant_starts.resize(nconstant);
    variable_starts.resize(nvariable);
    total_len=constant_lengths[0];

    for (size_t i=0; i<nvariable; ++i) {
        variable_starts[i]=total_len;
        total_len+=variable_lengths[i];
        constant_starts[i+1]=total_len;
        total_len+=constant_lengths[i+1];
    }

    return;
}

/* Searches a sequence and does some kind of operation to "OP". */

template<size_t nvariable, class OP>
bool search_sequence_internal(const char* seq, size_t len, 
    const std::vector<const sequence_dictionary*>& variable_dicts,
    const std::vector<std::u32string>& constant_hash,
    const std::vector<int>& variable_lengths,
    const std::vector<int>& constant_lengths,
    const std::vector<size_t>& variable_starts,
    const std::vector<size_t>& constant_starts,
    const size_t total_len,
    OP& operation) 
{
    if (len < total_len) {
        return false;
    } 

    constexpr size_t nconstant=nvariable+1;

    // Setting up the scanners for a perfect match or single mismatch.
    std::vector<hash_scanner> constant_scan;
    constant_scan.reserve(nconstant);
    for (size_t j=0; j<nconstant; ++j) {
        constant_scan.push_back(hash_scanner(seq+constant_starts[j], constant_lengths[j]));
    }

    std::vector<hash_scanner> variable_scan;
    variable_scan.reserve(nvariable);
    for (size_t j=0; j<nvariable; ++j) {
        variable_scan.push_back(hash_scanner(seq+variable_starts[j], variable_lengths[j]));
    }

    // Traversing through the sequence.
    size_t end=total_len;
    do {
        bool is_valid=true;
        for (auto& x : constant_scan) { is_valid &= x.valid(); } 
        for (auto& x : variable_scan) { is_valid &= x.valid(); }

        if (is_valid) {
            bool is_equal=true;
            for (size_t j=0; j<nconstant; ++j) {
                if (constant_hash[j]!=constant_scan[j].hash()) {
                    is_equal=false;
                    break;
                }
            }

            if (is_equal) {
                bool has_match=true;
                bool used_mismatch=false;
                combination<nvariable> tmp;

                for (size_t j=0; j<nvariable; ++j) {
                    const auto& curdict=*variable_dicts[j];
                    auto it=curdict.mapping.find(variable_scan[j].hash());

                    bool okay=false;
                    if (it!=curdict.mapping.end()) {
                        auto idx=(it->second).first;
                        if (idx!=-1) {
                            auto priority=(it->second).second;
                            if (priority==0 || !used_mismatch) {
                                if (priority==1) {
                                    used_mismatch=true; 
                                }
                                tmp.indices[j]=idx;
                                okay=true;
                            }
                        }
                    }

                    if (!okay) {
                        has_match=false;
                        break;
                    }
                }

                if (has_match) { 
                    operation(tmp);
                    return true;
                }
            }
        }

        if (end >= len) {
            return false;
        } 

        for (auto& x : constant_scan) { x.advance(); }
        for (auto& x : variable_scan) { x.advance(); }
        ++end;
    } while (1);
}

template<size_t nvariable, class OP>
bool search_sequence(const char* seq, size_t len, const search_info<nvariable>& info, OP& operation) 
{
    const auto& perfect=info.perfect;
    const auto& deleted=info.deleted;
    const auto& constant_hash=info.constant_hash;
    const auto& variable_lengths=info.variable_lengths;
    const auto& constant_lengths=info.constant_lengths;
    const auto& variable_starts=info.variable_starts;
    const auto& constant_starts=info.constant_starts;
    size_t total_len=info.total_len;
    bool allow_del=info.allow_del;

    std::vector<const sequence_dictionary*> ptrs(nvariable);
    for (size_t i=0; i<nvariable; ++i) { ptrs[i]=&(perfect[i]); }

    // Perfect matches or mismatch.
    if (search_sequence_internal<nvariable>(seq, len, 
        ptrs, constant_hash, 
        variable_lengths, constant_lengths,
        variable_starts, constant_starts,
        total_len, operation)) 
    {
        return true;
    } else if (!allow_del) {
        return false;
    }

    // Single deletion; introduced in one variable region at a time,
    // starting from the last.
    auto vlen=variable_lengths;
    const auto& clen=constant_lengths;
    auto vstart=variable_starts;
    auto cstart=constant_starts;

    for (size_t i=1; i<nvariable; ++i) { --vstart[i]; } // Preparing for deletion in the first variable region.
    for (size_t i=0; i<nvariable; ++i) { --cstart[i+1]; }

    for (size_t i=0; i<nvariable; ++i) {
        if (i) {
            ++vstart[i]; 
            ++cstart[i];
        }

        ptrs[i]=&(deleted[i]);
        --vlen[i];

        if (search_sequence_internal<nvariable>(seq, len, 
            ptrs, constant_hash, vlen, clen, vstart, cstart,
            total_len-1, operation))
        {
            return true;
        }

        ptrs[i]=&(perfect[i]);
        ++vlen[i];
    }

    return false;
}

#endif
