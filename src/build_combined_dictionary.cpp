#include "build_combined_dictionary.h"
#include "utils.h"

const std::vector<char> upper_bases = { 'A', 'C', 'G', 'T' };
const std::vector<char> lower_bases = { 'a', 'c', 'g', 't' };

void edit_sequence(combined_dictionary& ref, int n_sub, int n_insert, int n_del, int n_total, 
    const char* buffer, const size_t pos, const size_t len, 
    char* edited_buffer, size_t edited_pos, int nedits, 
    std::vector<int>& status)
{
    if (pos == len) {
        auto& mapping = ref[edited_pos];
        auto curhash = hash_sequence(edited_buffer, edited_pos);

        auto it = mapping.find(curhash);
        if (it != mapping.end()) {
            int& prev_nedits = (it->second).first;
            auto& prev_status = (it->second).second;

            if (prev_nedits == nedits) {
                if (prev_status != status) {
                    if (prev_nedits == 0) {
                        throw std::runtime_error("duplicated barcode sequences");
                    } else {
                        prev_nedits *= -1; // making it negative to indicate it shouldn't be trusted.
                    }
                }
            } else if (prev_nedits < 0) {
                if (-prev_nedits > nedits) {
                    prev_nedits = nedits;
                    prev_status = status;
                }
            } else {
                if (prev_nedits > nedits) {
                    prev_nedits = nedits;
                    prev_status = status;
                }
            }
        } else {
            mapping[curhash] = std::make_pair(nedits, status);
        }
        
        return;
    }

    const char curbuf = buffer[pos];
    edited_buffer[edited_pos] = curbuf;
    edit_sequence(ref, n_sub, n_insert, n_del, n_total,
        buffer, pos + 1, len,
        edited_buffer, edited_pos + 1, nedits,
        status);

    // Adding substitutions.
    if (n_sub && n_total) {
        for (size_t b = 0; b < upper_bases.size(); ++b) {
            if (curbuf != upper_bases[b] && curbuf != lower_bases[b]) {
                edited_buffer[edited_pos] = upper_bases[b];
                edit_sequence(ref, n_sub - 1, n_insert, n_del, n_total - 1,
                    buffer, pos + 1, len,
                    edited_buffer, edited_pos + 1, nedits + 1,
                    status);
            }
        }
        edited_buffer[edited_pos] = curbuf;
    }

    // Adding deletions.
    if (n_del && n_total) {
        edit_sequence(ref, n_sub, n_insert, n_del - 1, n_total - 1,
            buffer, pos + 1, len,
            edited_buffer, edited_pos, nedits + 1,
            status);
    }

    // Adding insertions. Don't worry about them cancelling out deletions, 
    // it'll get collapsed into the same entry in the unordered map anyway.
    //
    // Also note that the set-up means that we're adding insertions before 
    // the current position; in this case, we also ignore insertions at 
    // the start of the sequence, which don't matter to the final search. 
    if (n_insert && n_total && edited_pos!=0) {
        for (size_t b = 0; b < upper_bases.size(); ++b) {
            edited_buffer[edited_pos] = upper_bases[b];
            edit_sequence(ref, n_sub, n_insert - 1, n_del, n_total - 1,
                buffer, pos, len,
                edited_buffer, edited_pos + 1, nedits + 1,
                status);
        }
    }

    return;
}

void build_combined_dictionary_internal(combined_dictionary& ref, const Rcpp::StringVector& constant, const std::vector<Rcpp::StringVector>& variable, 
    int n_sub, int n_insert, int n_del, int n_total, 
    char* buffer, size_t pos, char* edit_buffer,
    std::vector<int>& status, int counter)
{
    // Fill in the preceding constant region (or the last constant region, if we're at the end).
    Rcpp::String con = constant[counter];
    const char* conptr = con.get_cstring();
    const size_t conlen = Rf_length(con.get_sexp());
    std::copy(conptr, conptr + conlen, buffer + pos);
    pos += conlen;

    if (counter == variable.size()) {
        // No need to add self, as edit_sequence will do it for us.
        edit_sequence(ref, n_sub, n_insert, n_del, n_total,
            buffer, 0, pos, // pos is now the length.
            edit_buffer, 0, 0, 
            status);
        return;
    }

    // Loop across all possibilities and fill them in.
    const auto& chosen = variable[counter];
    for (int i = 0; i < chosen.size(); ++i) {
        Rcpp::String var = chosen[i];
        const char* varptr = var.get_cstring();
        const size_t varlen = Rf_length(var.get_sexp());
        std::copy(varptr, varptr + varlen, buffer + pos);

        status[counter] = i;
        build_combined_dictionary_internal(ref, constant, variable, 
            n_sub, n_insert, n_del, n_total, 
            buffer, pos + varlen, edit_buffer, 
            status, counter + 1);
    }

    return;
}

combined_dictionary build_combined_dictionary(const Rcpp::StringVector& constant, const Rcpp::List& variable, int n_sub, int n_insert, int n_del, int n_total) { 
    std::vector<Rcpp::StringVector> varlist(variable.size());
    for (size_t v = 0; v < varlist.size(); ++v) { 
        varlist[v] = Rcpp::StringVector(variable[v]);
    }
    return build_combined_dictionary(constant, varlist, n_sub, n_insert, n_del, n_total);
}

combined_dictionary build_combined_dictionary(const Rcpp::StringVector& constant, const std::vector<Rcpp::StringVector>& variable, int n_sub, int n_insert, int n_del, int n_total) { 
    const size_t nvar = variable.size();
    if (constant.size() != nvar + 1) {
        throw std::runtime_error("number of constant regions must be 1 more than the number of variable regions");
    }

    if (n_sub < 0 || n_insert < 0 || n_del < 0 || n_total < 0) {
        throw std::runtime_error("maximum number of subsitutions, deletions, etc. must be a positive integer");
    }

    // Get full size here.
    size_t len = 0;
    for (size_t c = 0; c < constant.size(); ++c) {
        Rcpp::String con = constant[c];
        len += Rf_length(con.get_sexp());
    }

    for (auto& vars : variable) {
        len += check_length(vars);
    }

    // Instantiating the vectors.
    std::vector<char> buffer(len);
    std::vector<char> edit_buffer(len + n_insert);
    std::vector<int> status(nvar);

    combined_dictionary output;
    build_combined_dictionary_internal(output, constant, variable,
        n_sub, n_insert, n_del, n_total,
        buffer.data(), 0, edit_buffer.data(),
        status, 0);

    // Destroying all entries with ambiguous origins, as denoted by negative edits.
    for (auto& ep : output) {
        auto& e = ep.second;
        auto it = e.begin(); 
        while (it != e.end()) {
            if ((it->second).first < 0) {
                e.erase(it++);
            } else {
                ++it;
            }
        }
    }

    return output;
}
