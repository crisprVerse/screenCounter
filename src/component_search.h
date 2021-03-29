#ifndef COMPONENT_SEARCH_H 
#define COMPONENT_SEARCH_H

#include "Rcpp.h"
#include "build_basic_dictionary.h"
#include "hash_sequence.h"
#include "hash_scanner.h"
#include "utils.h"

/* Sets up search information. */

template<size_t nvariable>
struct combination {
    std::array<int, nvariable> indices;
};

template<size_t nvariable, size_t nconstant=nvariable+1>
struct search_info {
    search_info(Rcpp::List Guides, Rcpp::StringVector Constants) {   
        // Setting up the guides.
        if (nvariable!=Guides.size()) {
            std::stringstream err;
            err << "expecting " << nvariable << " variable regions";
            throw std::runtime_error(err.str());
        }

        for (size_t g=0; g<nvariable; ++g) {
            Rcpp::StringVector guides(Guides[g]);
            nbarcodes[g]=guides.size();

            size_t reflen = check_length(guides);
            variable_lengths[g]=reflen;

            perfect[g]=build_basic_dictionary(guides);
        }

        // Setting up the constant regions.
        if (Constants.size()!=nconstant) {
            throw std::runtime_error("number of constant regions should be 1 more than variable regions");
        }

        for (size_t s=0; s<nconstant; ++s) {
            Rcpp::String con(Constants[s]);
            const size_t len=Rf_length(con.get_sexp());
            const char* ptr=con.get_cstring();
            constant_hash[s]=hash_sequence(ptr, len); 
            constant_lengths[s]=len;
        }

        // Setting up search parameters.
        constant_starts[0]=0;
        total_len=constant_lengths[0];
        for (size_t i=0; i<nvariable; ++i) {
            variable_starts[i]=total_len;
            total_len+=variable_lengths[i];
            constant_starts[i+1]=total_len;
            total_len+=constant_lengths[i+1];
        }

        return;
    }

    std::array<basic_dictionary, nvariable> perfect;
    std::array<int, nvariable> variable_lengths;
    std::array<size_t, nvariable> variable_starts;
    std::array<size_t, nvariable> nbarcodes;

    std::array<seqhash, nconstant> constant_hash;
    std::array<int, nconstant> constant_lengths;
    std::array<size_t, nconstant> constant_starts;

    size_t total_len;
};

/* Searches a sequence and does some kind of operation to "OP". */

template<size_t nvariable, class OP, size_t nconstant=nvariable+1>
bool search_sequence_internal(const char* seq, size_t len, 
        
    const std::array<basic_dictionary, nvariable>& variable_dicts,
    const std::array<int, nvariable>& variable_lengths,
    const std::array<size_t, nvariable>& variable_starts,

    const std::array<seqhash, nconstant>& constant_hash,
    const std::array<int, nconstant>& constant_lengths,
    const std::array<size_t, nconstant>& constant_starts,

    const size_t total_len,
    OP& operation) 
{
    if (len < total_len) {
        return false;
    } 

    // Setting up the scanners for a perfect match or single mismatch.
    std::array<hash_scanner, nconstant> constant_scan;
    for (size_t j=0; j<nconstant; ++j) {
        constant_scan[j]=hash_scanner(seq+constant_starts[j], constant_lengths[j]);
    }

    std::array<hash_scanner, nvariable> variable_scan;
    for (size_t j=0; j<nvariable; ++j) {
        variable_scan[j]=hash_scanner(seq+variable_starts[j], variable_lengths[j]);
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
                combination<nvariable> tmp;

                for (size_t j=0; j<nvariable; ++j) {
                    const auto& curdict=variable_dicts[j];
                    auto it=curdict.find(variable_scan[j].hash());
                    if (it!=curdict.end()) {
                        tmp.indices[j]=it->second;
                    } else {
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

template<size_t nvariable, class OP, size_t nconstant=nvariable+1>
bool component_search(const char* seq, size_t len, const search_info<nvariable>& info, OP& operation) 
{
    const auto& perfect=info.perfect;
    const auto& constant_hash=info.constant_hash;
    const auto& variable_lengths=info.variable_lengths;
    const auto& constant_lengths=info.constant_lengths;
    const auto& variable_starts=info.variable_starts;
    const auto& constant_starts=info.constant_starts;
    size_t total_len=info.total_len;

    // Perfect matches or mismatch.
    if (search_sequence_internal<nvariable>(seq, len, 
        perfect, variable_lengths, variable_starts, 
        constant_hash, constant_lengths, constant_starts,
        total_len, operation)) 
    {
        return true;
    } else {
        return false;
    }
}

#endif
