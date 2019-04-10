#include "Rcpp.h"

#include "hash_sequence.h"
#include "build_hash.h"
#include <stdexcept>

/* Combination guide parser. */

// [[Rcpp::export(rng=false)]]
SEXP count_barcodes_combo(SEXP seqs, SEXP constants, SEXP guide_list) {
    // Setting up the guides.
    Rcpp::List Guides(guide_list);
    std::deque<seqhash> variable_hashes;
    std::deque<int> variable_lengths;

    for (size_t g=0; g<Guides.size(); ++g) {
        Rcpp::StringVector guides(Guides[g]);
        auto hash=build_hash(guides);
        variable_hashes.push_back(std::move(hash.first));
        
        if (hash.second.size()!=1) {
            throw std::runtime_error("all sequences should be of the same length");
        }
        variable_lengths.push_back(hash.second[0]);
    }

    // Setting up the constnt regions.
    Rcpp::StringVector Constants(constants);
    std::deque<std::u32string> constant_hash;
    std::deque<int> constant_lengths;

    for (size_t s=0; s<Constants.size(); ++s) {
        Rcpp::String con(Constants[s]);
        const size_t len=Rf_length(con.get_sexp());
        const char* ptr=con.get_cstring();
        constant_hash.push_back(hash_sequence(ptr, len)); 
        constant_lengths.push_back(len);
    }

    // Setting up search parameters.
    const size_t nvariable=variable_hashes.size(), 
        nconstant=constant_hash.size();
    if (nconstant!=nvariable+1) {
        throw std::runtime_error("number of constant regions should be 1 more than variable regions");
    }

    std::vector<size_t> constant_starts(nconstant), variable_starts(nvariable);
    size_t total_len=constant_lengths[0];
    for (size_t i=0; i<nvariable; ++i) {
        variable_starts[i]=total_len;
        total_len+=variable_lengths[i];
        constant_starts[i+1]=total_len;
        total_len+=constant_lengths[i+1];
    }

    // Running through the sequences and matching it to the guides.
    Rcpp::StringVector Seqs(seqs);
    std::vector<Rcpp::IntegerVector> output;
    for (size_t i=0; i<nvariable; ++i) { output[i]=Rcpp::IntegerVector(Seqs.size(), -1); }

    for (size_t i=0; i<Seqs.size(); ++i) {
        Rcpp::String s=Seqs[i];
        const char* sptr=s.get_cstring();
        const size_t len=Rf_length(s.get_sexp());
        if (len > total_len) { break; }

        // Setting up the scanners.
        std::vector<hash_scanner> constant_scan, variable_scan;
        constant_scan.reserve(nconstant);
        for (size_t i=0; i<nconstant; ++i) {
            constant_scan.push_back(hash_scanner(sptr+constant_starts[i], constant_lengths[i]));
        }
        variable_scan.reserve(nvariable);
        for (size_t i=0; i<nvariable; ++i) {
            variable_scan.push_back(hash_scanner(sptr+variable_starts[i], variable_lengths[i]));
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
                    if (constant_hash[i]!=constant_scan[i].hash()) {
                        is_equal=false;
                        break;
                    }
                }

                if (is_equal) {
                    for (size_t j=0; j<nvariable; ++j) {
                        const auto& curhash=variable_hashes[j];
                        auto it=curhash.find(variable_scan[j].hash());
                        if (it!=curhash.end()) {
                            output[j][i]=it->second;
                        }
                    }
                }
            }

            if (end >= len) {
                break;
            } 

            for (auto& x : constant_scan) { x.advance(); }
            for (auto& x : variable_scan) { x.advance(); }
            ++end;
        } while (1);
    }

    return Rcpp::List(output.begin(), output.end());
}
