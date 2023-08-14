#ifndef KAORI_PROCESS_DATA_HPP
#define KAORI_PROCESS_DATA_HPP

#include <thread>
#include "FastqReader.hpp"
#include "byteme/Reader.hpp"

/**
 * @file process_data.hpp
 *
 * @brief Process single- or paired-end data.
 */

namespace kaori {

/**
 * @cond
 */
template<bool use_names>
struct ChunkOfReads {
    ChunkOfReads() : sequence_offset(1), name_offset(1) {} // zero is always the first element.

    void clear() {
        sequence_buffer.clear();
        sequence_offset.resize(1);
        if constexpr(use_names) {
            name_buffer.clear();
            name_offset.resize(1);
        }
    }

    void add_read_sequence(const std::vector<char>& sequence) {
        add_read_details(sequence, sequence_buffer, sequence_offset);
    }

    void add_read_name(const std::vector<char>& name) {
        add_read_details(name, name_buffer, name_offset);
    }

    size_t size() const {
        return sequence_offset.size() - 1;
    }

    std::pair<const char*, const char*> get_sequence(size_t i) const {
        return get_details(i, sequence_buffer, sequence_offset);
    }

    std::pair<const char*, const char*> get_name(size_t i) const {
        return get_details(i, name_buffer, name_offset);
    }

private:
    std::vector<char> sequence_buffer;
    std::vector<size_t> sequence_offset;
    std::vector<char> name_buffer;
    std::vector<size_t> name_offset;

    static void add_read_details(const std::vector<char>& src, std::vector<char>& dst, std::vector<size_t>& offset) {
        dst.insert(dst.end(), src.begin(), src.end());
        auto last = offset.back();
        offset.push_back(last + src.size());
    }

    static std::pair<const char*, const char*> get_details(size_t i, const std::vector<char>& dest, const std::vector<size_t>& offset) {
        const char * base = dest.data();
        return std::make_pair(base + offset[i], base + offset[i + 1]);
    }
};
/**
 * @endcond
 */

/**
 * Run a handler for each read in single-end data.
 * This is done by calling `handler.process()` on each read.
 * It is expected that the results are stored in `handler` for retrieval by the caller.
 *
 * @tparam Handler A class that implements a handler for single-end data.
 *
 * @param input A `Reader` object containing data from a single-end FASTQ file.
 * @param handler Instance of the `Handler` class.
 * @param num_threads Number of threads to use for processing.
 * @param block_size Number of reads in each thread.
 *
 * @section single-handler-req Handler requirements
 * The `Handler` class is expected to implement the following methods:
 * - `initialize()`: this should be a `const` method that returns a state object (denoted here as having type `State`, though the exact name may vary).
 *   The idea is to store results in the state object for thread-safe execution.
 *   The state object should be default-constructible.
 * - `reduce(State& state)`: this should merge the results from the `state` object into the `Handler` instance.
 *   This will be called in a serial section and does not have to be thread-safe.
 *
 * The `Handler` should have a static `constexpr` variable `use_names`, indicating whether or not names should be passed to the `process()` method.
 *
 * If `use_names` is `false`, the `Handler` class should implement:
 * - `process(State& state, const std::pair<const char*, const char*>& seq)`: this should be a `const` method that processes the read in `seq` and stores its results in `state`.
 *   `seq` will contain pointers to the start and one-past-the-end of the read sequence.
 *
 * Otherwise, if `use_names` is `true`, the class should implement:
 * - `process(State& state, const std::pair<const char*, const char*>& name, const std::pair<const char*, const char*>& seq)`: 
 *    this should be a `const` method that processes the read in `seq` and stores its results in `state`.
 *   `name` will contain pointers to the start and one-past-the-end of the read name.
 *   `seq` will contain pointers to the start and one-past-the-end of the read sequence.
 */
template<class Handler>
void process_single_end_data(byteme::Reader* input, Handler& handler, int num_threads = 1, int block_size = 100000) {
    FastqReader fastq(input);
    bool finished = false;

    std::vector<ChunkOfReads<Handler::use_names> > reads(num_threads);
    std::vector<std::thread> jobs(num_threads);
    std::vector<decltype(handler.initialize())> states(num_threads);
    std::vector<std::string> errs(num_threads);

    auto join = [&](int i) -> void {
        if (jobs[i].joinable()) {
            jobs[i].join();
            if (errs[i] != "") {
                throw std::runtime_error(errs[i]);
            }
            handler.reduce(states[i]);
            reads[i].clear();
        }
    };

    // Safety measure to enforce const-ness within each thread.
    const Handler& conhandler = handler;

    try {
        while (!finished) {
            for (int t = 0; t < num_threads; ++t) {
                join(t);

                auto& curreads = reads[t];
                for (int b = 0; b < block_size; ++b) {
                    if (!fastq()) {
                        finished = true;
                        break;
                    }

                    curreads.add_read_sequence(fastq.get_sequence());
                    if constexpr(Handler::use_names) {
                        curreads.add_read_name(fastq.get_name());
                    }
                }

                states[t] = handler.initialize();
                jobs[t] = std::thread([&](int i) -> void {
                    try {
                        auto& state = states[i];
                        const auto& curreads = reads[i];
                        size_t nreads = curreads.size();

                        if constexpr(!Handler::use_names) {
                            for (size_t b = 0; b < nreads; ++b) {
                                conhandler.process(state, curreads.get_sequence(b));
                            }
                        } else {
                            for (size_t b = 0; b < nreads; ++b) {
                                conhandler.process(state, curreads.get_name(b), curreads.get_sequence(b));
                            }
                        }
                    } catch (std::exception& e) {
                        errs[i] = std::string(e.what());
                    }
                }, t);

                if (finished) {
                    // We won't get a future iteration to join the previous threads
                    // that we kicked off, so we do it now.
                    for (int u = 0; u < num_threads; ++u) {
                        auto pos = (u + t + 1) % num_threads;
                        join(pos);
                    }
                    break;
                }
            }
        }
    } catch (std::exception& e) {
        // Mopping up any loose threads, so to speak.
        for (int t = 0; t < num_threads; ++t) {
            if (jobs[t].joinable()) {
                jobs[t].join();
            }
        }
        throw;
    }

    return;
}

/**
 * Run a handler for each read in paired-end data, by calling `handler.process()` on each read pair.
 * It is expected that the results are stored in `handler` for retrieval by the caller.
 *
 * @tparam Handler A class that implements a handler for paired-end data.
 * @param input1 A `Reader` object containing data from the first FASTQ file in the pair.
 * @param input2 A `Reader` object containing data from the second FASTQ file in the pair.
 * @param handler Instance of the `Handler` class. 
 * @param num_threads Number of threads to use for processing.
 * @param block_size Number of reads in each thread.
 *
 * @section paired-handler-req Handler requirements
 * The `Handler` class is expected to implement the following methods:
 * - `initialize()`: this should be a `const` method that returns a state object (denoted here as having type `State`, though the exact name may vary).
 *   The idea is to store results in the state object for thread-safe execution.
 *   The state object should be default-constructible.
 * - `reduce(State& state)`: this should merge the results from the `state` object into the `Handler` instance.
 *   This will be called in a serial section and does not have to be thread-safe.
 *
 * The `Handler` should have a static `constexpr` variable `use_names`, indicating whether or not names should be passed to the `process()` method.
 *
 * If `use_names` is `false`, the `Handler` class should implement:
 * - `process(State& state, const std::pair<const char*, const char*>& seq1, const std::pair<const char*, const char*>& seq2)`: 
 *   this should be a `const` method that processes the paired reads in `seq1` and `seq2`, and stores its results in `state`.
 *   `seq1` and `seq2` will contain pointers to the start and one-past-the-end of the sequences of the paired reads.
 *
 * Otherwise, if `use_names` is `true`, the class should implement:
 * - `process(State& state, const std::pair<const char*, const char*>& name1, const std::pair<const char*, const char*>& seq1, const std::pair<const char*, const char*>& name2, const std::pair<const char*, const char*>& seq2)`: 
 *   this should be a `const` method that processes the paired reads in `seq1` and `seq2`, and stores its results in `state`.
 *   `name1` and `name2` will contain pointers to the start and one-past-the-end of the read names.
 *   `seq1` and `seq2` will contain pointers to the start and one-past-the-end of the read sequences.
 */
template<class Handler>
void process_paired_end_data(byteme::Reader* input1, byteme::Reader* input2, Handler& handler, int num_threads = 1, int block_size = 100000) {
    FastqReader fastq1(input1);
    FastqReader fastq2(input2);
    bool finished = false;

    std::vector<ChunkOfReads<Handler::use_names> > reads1(num_threads), reads2(num_threads);
    std::vector<std::thread> jobs(num_threads);
    std::vector<decltype(handler.initialize())> states(num_threads);
    std::vector<std::string> errs(num_threads);

    auto join = [&](int i) -> void {
        if (jobs[i].joinable()) {
            jobs[i].join();
            if (errs[i] != "") {
                throw std::runtime_error(errs[i]);
            }
            handler.reduce(states[i]);
            reads1[i].clear();
            reads2[i].clear();
        }
    };

    try {
        while (!finished) {
            for (int t = 0; t < num_threads; ++t) {
                join(t);

                bool finished1 = false;
                {
                    auto& curreads = reads1[t];
                    for (int b = 0; b < block_size; ++b) {
                        if (!fastq1()) {
                            finished1 = true;
                            break;
                        }

                        curreads.add_read_sequence(fastq1.get_sequence());
                        if constexpr(Handler::use_names) {
                            curreads.add_read_name(fastq1.get_name());
                        }
                    }
                }

                bool finished2 = false;
                {
                    auto& curreads = reads2[t];
                    for (int b = 0; b < block_size; ++b) {
                        if (!fastq2()) {
                            finished2 = true;
                            break;
                        }

                        curreads.add_read_sequence(fastq2.get_sequence());
                        if constexpr(Handler::use_names) {
                            curreads.add_read_name(fastq2.get_name());
                        }
                    }
                }

                if (finished1 != finished2 || reads1[t].size() != reads2[t].size()) {
                    throw std::runtime_error("different number of reads in paired FASTQ files");
                } else if (finished1) {
                    finished = true;
                }

                states[t] = handler.initialize();
                jobs[t] = std::thread([&](int i) -> void {
                    auto& state = states[i];
                    const auto& curreads1 = reads1[i];
                    const auto& curreads2 = reads2[i];
                    size_t nreads = curreads1.size();

                    try {
                        if constexpr(!Handler::use_names) {
                            for (size_t b = 0; b < nreads; ++b) {
                                handler.process(state, curreads1.get_sequence(b), curreads2.get_sequence(b));
                            }
                        } else {
                            for (size_t b = 0; b < nreads; ++b) {
                                handler.process(
                                    state,
                                    curreads1.get_name(b), 
                                    curreads1.get_sequence(b),
                                    curreads2.get_name(b), 
                                    curreads2.get_sequence(b)
                                );
                            }
                        }
                    } catch (std::exception& e) {
                        errs[i] = std::string(e.what());
                    }
                }, t);

                if (finished) {
                    // We won't get a future iteration to join the previous threads
                    // that we kicked off, so we do it now.
                    for (int u = 0; u < num_threads; ++u) {
                        auto pos = (u + t + 1) % num_threads;
                        join(pos);
                    }
                    break;
                }
            }
        }
    } catch (std::exception& e) {
        // Mopping up any loose threads, so to speak.
        for (int t = 0; t < num_threads; ++t) {
            if (jobs[t].joinable()) {
                jobs[t].join();
            }
        }
        throw;
    }

    return;
}


}

#endif
