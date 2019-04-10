#include "Rcpp.h"

#include "hash_sequence.h"
#include "build_hash.h"
#include <stdexcept>
#include <vector>

/* Combination guide parser. */

struct se_combo_info {
    se_combo_info(Rcpp::List, Rcpp::StringVector);
    std::vector<seqhash> variable_hashes;
    std::vector<std::u32string> constant_hash;
    std::vector<int> constant_lengths, variable_lengths;
    std::vector<size_t> constant_starts, variable_starts;
    size_t total_len;
};

se_combo_info::se_combo_info(Rcpp::List Guides, Rcpp::StringVector Constants) {
    // Setting up the guides.
    size_t nvariable=Guides.size();
    variable_hashes.reserve(nvariable);
    variable_lengths.reserve(nvariable);

    for (size_t g=0; g<nvariable; ++g) {
        Rcpp::StringVector guides(Guides[g]);
        auto hash=build_hash(guides);
        variable_hashes.push_back(std::move(hash.first));
        
        if (hash.second.size()!=1) {
            throw std::runtime_error("all sequences should be of the same length");
        }
        variable_lengths.push_back(hash.second[0]);
    }

    // Setting up the constnt regions.
    size_t nconstant=Constants.size();
    if (nconstant!=nvariable+1) {
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
    constant_starts.reserve(nconstant);
    variable_starts.reserve(nvariable);
    total_len=constant_lengths[0];

    for (size_t i=0; i<nvariable; ++i) {
        variable_starts[i]=total_len;
        total_len+=variable_lengths[i];
        constant_starts[i+1]=total_len;
        total_len+=constant_lengths[i+1];
    }

    return;
}

// [[Rcpp::export(rng=false)]]
SEXP setup_barcodes_combo(SEXP constants, SEXP guide_list) {
    se_combo_info * ptr = new se_combo_info(guide_list, constants);
    return Rcpp::XPtr<se_combo_info>(ptr, true);
}

// [[Rcpp::export(rng=false)]]
SEXP count_barcodes_combo(SEXP seqs, SEXP xptr) {
    Rcpp::XPtr<se_combo_info> ptr(xptr);
    const auto& variable_hashes=ptr->variable_hashes;
    const auto& constant_hash=ptr->constant_hash;
    const auto& variable_lengths=ptr->variable_lengths;
    const auto& constant_lengths=ptr->constant_lengths;
    const auto& variable_starts=ptr->variable_starts;
    const auto& constant_starts=ptr->constant_starts;

    const size_t total_len=ptr->total_len;
    const size_t nconstant=constant_hash.size();
    const size_t nvariable=variable_hashes.size();

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
