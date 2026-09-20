// Exact target-weight histogram for ShiftAlpha64 applied to all EBCH weight-22 messages.
//
// The hot loop advances both the GF(2^64) message and its 128-bit extended-BCH
// codeword.  For the fixed systematic encoder B and T(m)=gamma*m,
//
//   B(Tm) = shift(B(m) + m_63 B(e_63))
//           + parity(m mod 2^63)e_127 + m_63 B(0x1b).
//
// This identity avoids eight table lookups per source and exponent.  Threads
// own disjoint contiguous source ranges for the complete scan.

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

constexpr std::uint64_t kFieldReduction = 0x1bULL;
constexpr std::uint64_t kLower63Mask = 0x7fffffffffffffffULL;
constexpr std::uint64_t kGenerator63Low = 0x8000000000000000ULL;
constexpr std::uint64_t kGenerator63High = 0xfa422a8c5cac150fULL;
constexpr std::uint64_t kReductionCodeLow = 0xf1eb066794ab8f29ULL;
constexpr std::uint64_t kReductionCodeHigh = 0x0000000000000009ULL;

struct Codeword {
    std::uint64_t low;
    std::uint64_t high;
};

constexpr std::array<Codeword, 64> kGeneratorRows{{
    {0xf4845518b9582a1fULL,0x8000000000000000ULL},
    {0xe908aa3172b0543eULL,0x8000000000000001ULL},
    {0xd2115462e560a87cULL,0x8000000000000003ULL},
    {0xa422a8c5cac150f8ULL,0x8000000000000007ULL},
    {0x4845518b9582a1f0ULL,0x800000000000000fULL},
    {0x908aa3172b0543e0ULL,0x800000000000001eULL},
    {0x2115462e560a87c0ULL,0x800000000000003dULL},
    {0x422a8c5cac150f80ULL,0x800000000000007aULL},
    {0x845518b9582a1f00ULL,0x80000000000000f4ULL},
    {0x08aa3172b0543e00ULL,0x80000000000001e9ULL},
    {0x115462e560a87c00ULL,0x80000000000003d2ULL},
    {0x22a8c5cac150f800ULL,0x80000000000007a4ULL},
    {0x45518b9582a1f000ULL,0x8000000000000f48ULL},
    {0x8aa3172b0543e000ULL,0x8000000000001e90ULL},
    {0x15462e560a87c000ULL,0x8000000000003d21ULL},
    {0x2a8c5cac150f8000ULL,0x8000000000007a42ULL},
    {0x5518b9582a1f0000ULL,0x800000000000f484ULL},
    {0xaa3172b0543e0000ULL,0x800000000001e908ULL},
    {0x5462e560a87c0000ULL,0x800000000003d211ULL},
    {0xa8c5cac150f80000ULL,0x800000000007a422ULL},
    {0x518b9582a1f00000ULL,0x80000000000f4845ULL},
    {0xa3172b0543e00000ULL,0x80000000001e908aULL},
    {0x462e560a87c00000ULL,0x80000000003d2115ULL},
    {0x8c5cac150f800000ULL,0x80000000007a422aULL},
    {0x18b9582a1f000000ULL,0x8000000000f48455ULL},
    {0x3172b0543e000000ULL,0x8000000001e908aaULL},
    {0x62e560a87c000000ULL,0x8000000003d21154ULL},
    {0xc5cac150f8000000ULL,0x8000000007a422a8ULL},
    {0x8b9582a1f0000000ULL,0x800000000f484551ULL},
    {0x172b0543e0000000ULL,0x800000001e908aa3ULL},
    {0x2e560a87c0000000ULL,0x800000003d211546ULL},
    {0x5cac150f80000000ULL,0x800000007a422a8cULL},
    {0xb9582a1f00000000ULL,0x80000000f4845518ULL},
    {0x72b0543e00000000ULL,0x80000001e908aa31ULL},
    {0xe560a87c00000000ULL,0x80000003d2115462ULL},
    {0xcac150f800000000ULL,0x80000007a422a8c5ULL},
    {0x9582a1f000000000ULL,0x8000000f4845518bULL},
    {0x2b0543e000000000ULL,0x8000001e908aa317ULL},
    {0x560a87c000000000ULL,0x8000003d2115462eULL},
    {0xac150f8000000000ULL,0x8000007a422a8c5cULL},
    {0x582a1f0000000000ULL,0x800000f4845518b9ULL},
    {0xb0543e0000000000ULL,0x800001e908aa3172ULL},
    {0x60a87c0000000000ULL,0x800003d2115462e5ULL},
    {0xc150f80000000000ULL,0x800007a422a8c5caULL},
    {0x82a1f00000000000ULL,0x80000f4845518b95ULL},
    {0x0543e00000000000ULL,0x80001e908aa3172bULL},
    {0x0a87c00000000000ULL,0x80003d2115462e56ULL},
    {0x150f800000000000ULL,0x80007a422a8c5cacULL},
    {0x2a1f000000000000ULL,0x8000f4845518b958ULL},
    {0x543e000000000000ULL,0x8001e908aa3172b0ULL},
    {0xa87c000000000000ULL,0x8003d2115462e560ULL},
    {0x50f8000000000000ULL,0x8007a422a8c5cac1ULL},
    {0xa1f0000000000000ULL,0x800f4845518b9582ULL},
    {0x43e0000000000000ULL,0x801e908aa3172b05ULL},
    {0x87c0000000000000ULL,0x803d2115462e560aULL},
    {0x0f80000000000000ULL,0x807a422a8c5cac15ULL},
    {0x1f00000000000000ULL,0x80f4845518b9582aULL},
    {0x3e00000000000000ULL,0x81e908aa3172b054ULL},
    {0x7c00000000000000ULL,0x83d2115462e560a8ULL},
    {0xf800000000000000ULL,0x87a422a8c5cac150ULL},
    {0xf000000000000000ULL,0x8f4845518b9582a1ULL},
    {0xe000000000000000ULL,0x9e908aa3172b0543ULL},
    {0xc000000000000000ULL,0xbd2115462e560a87ULL},
    {0x8000000000000000ULL,0xfa422a8c5cac150fULL},
}};

Codeword encode(std::uint64_t message) {
    Codeword result{0, 0};
    while (message != 0) {
        const unsigned bit = std::countr_zero(message);
        result.low ^= kGeneratorRows[bit].low;
        result.high ^= kGeneratorRows[bit].high;
        message &= message - 1;
    }
    return result;
}

unsigned weight(const Codeword codeword) {
    return std::popcount(codeword.low) + std::popcount(codeword.high);
}

void advance(std::uint64_t& message, Codeword& codeword) {
    const std::uint64_t top = message >> 63;
    const std::uint64_t mask = 0ULL - top;
    const std::uint64_t lower = message & kLower63Mask;
    const std::uint64_t old_low = codeword.low ^ (kGenerator63Low & mask);
    const std::uint64_t old_high = codeword.high ^ (kGenerator63High & mask);
    codeword.low = old_low << 1;
    codeword.high = (old_high << 1) | (old_low >> 63);
    codeword.high ^= static_cast<std::uint64_t>(std::popcount(lower) & 1U) << 63;
    codeword.low ^= kReductionCodeLow & mask;
    codeword.high ^= kReductionCodeHigh & mask;
    message = (message << 1) ^ (kFieldReduction & mask);
}

std::vector<std::uint64_t> read_messages(const std::filesystem::path& path) {
    const auto bytes = std::filesystem::file_size(path);
    if (bytes == 0 || bytes % sizeof(std::uint64_t) != 0) {
        throw std::runtime_error("source-message file has invalid size");
    }
    std::vector<std::uint64_t> messages(
        static_cast<std::size_t>(bytes / sizeof(std::uint64_t)));
    std::ifstream input(path, std::ios::binary);
    input.read(reinterpret_cast<char*>(messages.data()), static_cast<std::streamsize>(bytes));
    if (!input || input.gcount() != static_cast<std::streamsize>(bytes)) {
        throw std::runtime_error("failed to read the complete minimum-message file");
    }
    if (!std::is_sorted(messages.begin(), messages.end()) ||
        std::adjacent_find(messages.begin(), messages.end()) != messages.end() ||
        messages.front() == 0) {
        throw std::runtime_error("source messages are not sorted unique nonzero values");
    }
    return messages;
}

struct WorkerResult {
    std::array<std::uint64_t, 129> histogram{};
    std::vector<std::uint64_t> low_by_exponent;
    unsigned minimum_weight = 129;
};

std::uint64_t spectrum_count(unsigned target_weight) {
    switch (target_weight) {
        case 22: return 243840ULL;
        case 24: return 6855968ULL;
        case 26: return 107988608ULL;
        case 28: return 1479751168ULL;
        case 30: return 16581217536ULL;
        case 32: return 161471882796ULL;
        default: return 0;
    }
}

}  // namespace

int main(int argc, char** argv) try {
    if (argc < 2 || argc > 7) {
        std::cerr << "usage: scanner messages.bin [lag_start=64] [lag_count=16384] "
                     "[threads=hardware] [low_report_max=30] [source_weight=22]\n";
        return 2;
    }
    const std::filesystem::path path = argv[1];
    const unsigned lag_start = argc >= 3 ? static_cast<unsigned>(std::stoul(argv[2])) : 64U;
    const unsigned lag_count = argc >= 4 ? static_cast<unsigned>(std::stoul(argv[3])) : 16384U;
    const unsigned hardware = std::max(1U, std::thread::hardware_concurrency());
    const unsigned requested_threads = argc >= 5 ? static_cast<unsigned>(std::stoul(argv[4])) : hardware;
    const unsigned thread_count = std::max(1U, std::min(requested_threads, 16U));
    const unsigned low_report_max = argc >= 6 ? static_cast<unsigned>(std::stoul(argv[5])) : 30U;
    const unsigned source_weight = argc >= 7 ? static_cast<unsigned>(std::stoul(argv[6])) : 22U;
    if (lag_count == 0 || low_report_max > 128) {
        throw std::runtime_error("lag count must be positive and low report maximum at most 128");
    }

    auto messages = read_messages(path);
    std::vector<Codeword> codewords(messages.size());
    for (std::size_t index = 0; index < messages.size(); ++index) {
        codewords[index] = encode(messages[index]);
        if (weight(codewords[index]) != source_weight) {
            throw std::runtime_error("input contains a BCH message outside the declared source shell");
        }
    }
    for (unsigned step = 0; step < lag_start; ++step) {
        for (std::size_t index = 0; index < messages.size(); ++index) {
            advance(messages[index], codewords[index]);
        }
    }
    for (std::size_t probe = 0; probe < 17; ++probe) {
        const std::size_t index = probe * (messages.size() - 1) / 16;
        const Codeword direct = encode(messages[index]);
        if (direct.low != codewords[index].low || direct.high != codewords[index].high) {
            throw std::runtime_error("rolling BCH recurrence failed direct re-encoding validation");
        }
    }

    std::vector<WorkerResult> results(thread_count);
    for (auto& result : results) {
        result.low_by_exponent.assign(
            static_cast<std::size_t>(lag_count) * (low_report_max + 1), 0);
    }
    std::atomic<unsigned> progress{0};
    const auto started = std::chrono::steady_clock::now();
    std::vector<std::thread> workers;
    workers.reserve(thread_count);
    for (unsigned thread_index = 0; thread_index < thread_count; ++thread_index) {
        const std::size_t begin = messages.size() * thread_index / thread_count;
        const std::size_t end = messages.size() * (thread_index + 1) / thread_count;
        workers.emplace_back([&, thread_index, begin, end]() {
            WorkerResult& result = results[thread_index];
            for (unsigned lag_offset = 0; lag_offset < lag_count; ++lag_offset) {
                const std::size_t low_base = static_cast<std::size_t>(lag_offset) * (low_report_max + 1);
                for (std::size_t index = begin; index < end; ++index) {
                    const unsigned target_weight = weight(codewords[index]);
                    ++result.histogram[target_weight];
                    result.minimum_weight = std::min(result.minimum_weight, target_weight);
                    if (target_weight <= low_report_max) {
                        ++result.low_by_exponent[low_base + target_weight];
                    }
                    advance(messages[index], codewords[index]);
                }
                if (thread_index == 0 && ((lag_offset + 1) % 1024 == 0)) {
                    progress.store(lag_offset + 1, std::memory_order_relaxed);
                    const double elapsed = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - started).count();
                    std::cerr << "scanned " << (lag_offset + 1) << '/' << lag_count
                              << " exponents in " << std::fixed << std::setprecision(1)
                              << elapsed << "s\n";
                }
            }
        });
    }
    for (auto& worker : workers) {
        worker.join();
    }
    for (std::size_t probe = 0; probe < 17; ++probe) {
        const std::size_t index = probe * (messages.size() - 1) / 16;
        const Codeword direct = encode(messages[index]);
        if (direct.low != codewords[index].low || direct.high != codewords[index].high) {
            throw std::runtime_error("rolling BCH recurrence failed final direct re-encoding validation");
        }
    }
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();

    std::array<std::uint64_t, 129> histogram{};
    unsigned minimum_weight = 129;
    std::vector<std::uint64_t> low_by_exponent(
        static_cast<std::size_t>(lag_count) * (low_report_max + 1), 0);
    for (const auto& result : results) {
        minimum_weight = std::min(minimum_weight, result.minimum_weight);
        for (unsigned target_weight = 0; target_weight <= 128; ++target_weight) {
            histogram[target_weight] += result.histogram[target_weight];
        }
        for (std::size_t index = 0; index < low_by_exponent.size(); ++index) {
            low_by_exponent[index] += result.low_by_exponent[index];
        }
    }
    const std::uint64_t expected_pairs =
        static_cast<std::uint64_t>(messages.size()) * lag_count;
    std::uint64_t observed_pairs = 0;
    for (const auto count : histogram) {
        observed_pairs += count;
    }
    if (observed_pairs != expected_pairs) {
        throw std::runtime_error("target histogram failed total-mass validation");
    }
    for (unsigned odd_weight = 1; odd_weight <= 127; odd_weight += 2) {
        if (histogram[odd_weight] != 0) {
            throw std::runtime_error("extended BCH histogram contains an odd weight");
        }
    }

    std::cout << std::setprecision(17);
    std::cout << "{\n"
              << "  \"schema\": \"riffle-shiftalpha64-source-shell-target-spectrum-v1\",\n"
              << "  \"evidence_label\": \"EXACT_COMPLETE_SOURCE_SHELL_TARGET_HISTOGRAM\",\n"
              << "  \"construction\": \"Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4\",\n"
              << "  \"source_bch_weight\": " << source_weight << ",\n"
              << "  \"source_count\": " << messages.size() << ",\n"
              << "  \"exponent_start_inclusive\": " << lag_start << ",\n"
              << "  \"exponent_end_exclusive\": " << (lag_start + lag_count) << ",\n"
              << "  \"source_exponent_pairs\": " << expected_pairs << ",\n"
              << "  \"minimum_target_bch_weight\": " << minimum_weight << ",\n"
              << "  \"minimum_one_symbol_profile_total_weight\": "
              << (2 * source_weight + minimum_weight) << ",\n"
              << "  \"target_weight_histogram\": {\n";
    bool first = true;
    for (unsigned target_weight = 0; target_weight <= 128; ++target_weight) {
        if (histogram[target_weight] == 0) {
            continue;
        }
        if (!first) {
            std::cout << ",\n";
        }
        first = false;
        std::cout << "    \"" << target_weight << "\": " << histogram[target_weight];
    }
    std::cout << "\n  },\n  \"low_weight_hits_by_exponent\": [";
    bool first_row = true;
    for (unsigned lag_offset = 0; lag_offset < lag_count; ++lag_offset) {
        for (unsigned target_weight = 0; target_weight <= low_report_max; ++target_weight) {
            const auto count = low_by_exponent[
                static_cast<std::size_t>(lag_offset) * (low_report_max + 1) + target_weight];
            if (count == 0) {
                continue;
            }
            std::cout << (first_row ? "\n" : ",\n")
                      << "    {\"exponent\": " << (lag_start + lag_offset)
                      << ", \"target_weight\": " << target_weight
                      << ", \"ordered_pairs\": " << count << "}";
            first_row = false;
        }
    }
    if (!first_row) {
        std::cout << '\n';
    }
    std::cout << "  ],\n  \"uniform_alpha_comparison\": {\n";
    first = true;
    const long double field_order = 18446744073709551615.0L;
    for (unsigned target_weight : {22U, 24U, 26U, 28U, 30U, 32U}) {
        const long double expected = static_cast<long double>(expected_pairs)
            * static_cast<long double>(spectrum_count(target_weight)) / field_order;
        if (!first) {
            std::cout << ",\n";
        }
        first = false;
        std::cout << "    \"" << target_weight << "\": {\"expected_pairs\": "
                  << static_cast<double>(expected) << ", \"observed_pairs\": "
                  << histogram[target_weight] << "}";
    }
    std::cout << "\n  },\n"
              << "  \"validation\": {\n"
              << "    \"input_sorted_unique_nonzero\": true,\n"
              << "    \"all_source_bch_weights_equal_declared_weight\": true,\n"
              << "    \"rolling_recurrence_matches_direct_reencoding_at_both_endpoints\": true,\n"
              << "    \"histogram_total_mass_passes\": true,\n"
              << "    \"all_target_weights_even\": true\n"
              << "  },\n"
              << "  \"threads\": " << thread_count << ",\n"
              << "  \"elapsed_seconds\": " << elapsed << ",\n"
              << "  \"scope\": \"Exact only for the declared source BCH weight.\"\n"
              << "}\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
}
