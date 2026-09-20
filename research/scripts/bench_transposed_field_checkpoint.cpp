#define main bench_transposed_structured_step_hidden_main
#include "bench_transposed_structured_step.cpp"
#undef main

namespace {

enum class CheckpointImplementation {
    accumulator,
    direct,
    bitsliced,
    rectangular,
    rectangular_avx2,
};

constexpr std::size_t checkpoint_epoch_blocks = 256;
constexpr unsigned checkpoint_state_blocks = 64;

inline void transpose64_standard(std::uint64_t* rows) noexcept {
    std::uint64_t mask = 0x00000000ffffffffULL;
    for (unsigned shift = 32; shift != 0; shift >>= 1) {
        for (unsigned base = 0; base < 64; base = (base + shift + 1) & ~shift) {
            const std::uint64_t difference =
                ((rows[base] >> shift) ^ rows[base + shift]) & mask;
            rows[base] ^= difference << shift;
            rows[base + shift] ^= difference;
        }
        mask ^= mask << (shift >> 1);
    }
}

template <int Shift, std::uint64_t Mask>
inline void transpose64_avx2_large_stage(std::uint64_t* rows) noexcept {
    const __m256i mask = _mm256_set1_epi64x(static_cast<long long>(Mask));
    for (unsigned base = 0; base < 64; base += 2 * Shift) {
        for (unsigned offset = 0; offset < Shift; offset += 4) {
            __m256i low = _mm256_load_si256(
                reinterpret_cast<const __m256i*>(rows + base + offset));
            __m256i high = _mm256_load_si256(
                reinterpret_cast<const __m256i*>(rows + base + Shift + offset));
            const __m256i difference = _mm256_and_si256(
                _mm256_xor_si256(_mm256_srli_epi64(low, Shift), high), mask);
            low = _mm256_xor_si256(low, _mm256_slli_epi64(difference, Shift));
            high = _mm256_xor_si256(high, difference);
            _mm256_store_si256(
                reinterpret_cast<__m256i*>(rows + base + offset), low);
            _mm256_store_si256(
                reinterpret_cast<__m256i*>(rows + base + Shift + offset), high);
        }
    }
}

inline void transpose64_standard_avx2(std::uint64_t* rows) noexcept {
    transpose64_avx2_large_stage<32, 0x00000000ffffffffULL>(rows);
    transpose64_avx2_large_stage<16, 0x0000ffff0000ffffULL>(rows);
    transpose64_avx2_large_stage<8, 0x00ff00ff00ff00ffULL>(rows);
    transpose64_avx2_large_stage<4, 0x0f0f0f0f0f0f0f0fULL>(rows);

    const __m256i mask2 = _mm256_set1_epi64x(0x3333333333333333ULL);
    for (unsigned base = 0; base < 64; base += 4) {
        __m256i rows4 = _mm256_load_si256(
            reinterpret_cast<const __m256i*>(rows + base));
        const __m256i partner = _mm256_permute4x64_epi64(rows4, 0x4e);
        __m256i difference = _mm256_and_si256(
            _mm256_xor_si256(_mm256_srli_epi64(rows4, 2), partner), mask2);
        difference = _mm256_permute4x64_epi64(difference, 0x44);
        const __m256i adjustment = _mm256_blend_epi32(
            _mm256_slli_epi64(difference, 2), difference, 0xf0);
        rows4 = _mm256_xor_si256(rows4, adjustment);
        _mm256_store_si256(reinterpret_cast<__m256i*>(rows + base), rows4);
    }

    const __m256i mask1 = _mm256_set1_epi64x(0x5555555555555555ULL);
    for (unsigned base = 0; base < 64; base += 4) {
        __m256i rows4 = _mm256_load_si256(
            reinterpret_cast<const __m256i*>(rows + base));
        const __m256i partner = _mm256_permute4x64_epi64(rows4, 0xb1);
        __m256i difference = _mm256_and_si256(
            _mm256_xor_si256(_mm256_srli_epi64(rows4, 1), partner), mask1);
        difference = _mm256_permute4x64_epi64(difference, 0xa0);
        const __m256i adjustment = _mm256_blend_epi32(
            _mm256_slli_epi64(difference, 1), difference, 0xcc);
        rows4 = _mm256_xor_si256(rows4, adjustment);
        _mm256_store_si256(reinterpret_cast<__m256i*>(rows + base), rows4);
    }
}

template <void (*Transpose)(std::uint64_t*)>
inline void transform_field_rectangle64(
    Block* state, std::uint64_t coefficient) noexcept {
    alignas(64) std::array<std::uint64_t, 128> lanes;
    for (unsigned row = 0; row < checkpoint_state_blocks; ++row) {
        lanes[row] = static_cast<std::uint64_t>(
            _mm_cvtsi128_si64(state[row].value));
        lanes[64 + row] = static_cast<std::uint64_t>(
            _mm_extract_epi64(state[row].value, 1));
    }
    Transpose(lanes.data());
    Transpose(lanes.data() + 64);

    for (unsigned bit = 0; bit < 64; ++bit) {
        const __m256i packed = _mm256_set_epi64x(
            0,
            static_cast<long long>(lanes[64 + bit]),
            0,
            static_cast<long long>(lanes[bit]));
        const __m256i product =
            field_multiply_vector<checkpoint_state_blocks>(packed, coefficient);
        lanes[bit] = static_cast<std::uint64_t>(
            _mm_cvtsi128_si64(_mm256_castsi256_si128(product)));
        lanes[64 + bit] = static_cast<std::uint64_t>(
            _mm_cvtsi128_si64(_mm256_extracti128_si256(product, 1)));
    }

    Transpose(lanes.data());
    Transpose(lanes.data() + 64);
    for (unsigned row = 0; row < checkpoint_state_blocks; ++row) {
        state[row].value = _mm_set_epi64x(
            static_cast<long long>(lanes[64 + row]),
            static_cast<long long>(lanes[row]));
    }
}

template <CheckpointImplementation Implementation>
void apply_checkpoint_chain(
    Block* word,
    std::size_t word_blocks,
    const std::uint64_t* coefficients) noexcept {
    std::array<Block, checkpoint_state_blocks> state;
    state.fill(zero_block());
    alignas(64) std::array<Block, 128> square{};
    Block* active_state = state.data();
    if constexpr (Implementation == CheckpointImplementation::bitsliced) {
        active_state = square.data();
    }
    const std::size_t epochs = word_blocks / checkpoint_epoch_blocks;

    for (std::size_t epoch = epochs; epoch-- > 0;) {
        if (epoch + 1 < epochs) {
            if constexpr (Implementation == CheckpointImplementation::direct) {
                apply_direct_field<checkpoint_state_blocks>(
                    active_state, coefficients[epoch]);
            } else if constexpr (
                Implementation == CheckpointImplementation::bitsliced) {
                transform_field_square<checkpoint_state_blocks>(
                    square.data(), coefficients[epoch]);
            } else if constexpr (
                Implementation == CheckpointImplementation::rectangular) {
                transform_field_rectangle64<transpose64_standard>(
                    active_state, coefficients[epoch]);
            } else if constexpr (
                Implementation == CheckpointImplementation::rectangular_avx2) {
                transform_field_rectangle64<transpose64_standard_avx2>(
                    active_state, coefficients[epoch]);
            }
        }

        Block* const position = word + epoch * checkpoint_epoch_blocks;
        for (std::size_t offset = checkpoint_epoch_blocks; offset-- > 0;) {
            const unsigned lane = static_cast<unsigned>(
                offset & (checkpoint_state_blocks - 1));
            active_state[lane] = xor_block(active_state[lane], position[offset]);
            position[offset] = active_state[lane];
        }
    }
}

const char* implementation_name(CheckpointImplementation implementation) {
    switch (implementation) {
        case CheckpointImplementation::accumulator:
            return "accumulator_only";
        case CheckpointImplementation::direct:
            return "block_four_russians_compact_coefficients";
        case CheckpointImplementation::bitsliced:
            return "avx2_transpose128_vpclmulqdq_two_lanes";
        case CheckpointImplementation::rectangular:
            return "transpose64x128_vpclmulqdq_two_lanes";
        case CheckpointImplementation::rectangular_avx2:
            return "avx2_transpose64x128_vpclmulqdq_two_lanes";
    }
    return "unknown";
}

CheckpointImplementation parse_checkpoint_implementation(const std::string& name) {
    if (name == "accumulator") {
        return CheckpointImplementation::accumulator;
    }
    if (name == "direct") {
        return CheckpointImplementation::direct;
    }
    if (name == "bitsliced") {
        return CheckpointImplementation::bitsliced;
    }
    if (name == "rectangular") {
        return CheckpointImplementation::rectangular;
    }
    if (name == "rectangular-avx2") {
        return CheckpointImplementation::rectangular_avx2;
    }
    throw std::invalid_argument(
        "implementation must be accumulator, direct, bitsliced, rectangular, or rectangular-avx2");
}

void invoke_checkpoint(
    CheckpointImplementation implementation,
    Block* word,
    std::size_t word_blocks,
    const std::uint64_t* coefficients) noexcept {
    switch (implementation) {
        case CheckpointImplementation::accumulator:
            apply_checkpoint_chain<CheckpointImplementation::accumulator>(
                word, word_blocks, coefficients);
            return;
        case CheckpointImplementation::direct:
            apply_checkpoint_chain<CheckpointImplementation::direct>(
                word, word_blocks, coefficients);
            return;
        case CheckpointImplementation::bitsliced:
            apply_checkpoint_chain<CheckpointImplementation::bitsliced>(
                word, word_blocks, coefficients);
            return;
        case CheckpointImplementation::rectangular:
            apply_checkpoint_chain<CheckpointImplementation::rectangular>(
                word, word_blocks, coefficients);
            return;
        case CheckpointImplementation::rectangular_avx2:
            apply_checkpoint_chain<CheckpointImplementation::rectangular_avx2>(
                word, word_blocks, coefficients);
            return;
    }
}

void self_test_checkpoint() {
    std::uint64_t transpose_random_state = 0x7472616e73706f73ULL;
    for (unsigned trial = 0; trial < 256; ++trial) {
        alignas(64) std::array<std::uint64_t, 64> scalar;
        alignas(64) std::array<std::uint64_t, 64> vector;
        for (std::uint64_t& value : scalar) {
            value = splitmix64(transpose_random_state);
        }
        vector = scalar;
        const auto original = scalar;
        transpose64_standard(scalar.data());
        transpose64_standard_avx2(vector.data());
        if (scalar != vector) {
            throw std::runtime_error(
                "AVX2 64-by-64 transpose disagrees with scalar reference");
        }
        transpose64_standard_avx2(vector.data());
        if (vector != original) {
            throw std::runtime_error("AVX2 64-by-64 transpose is not involutive");
        }
    }

    constexpr std::size_t blocks = 3 * checkpoint_epoch_blocks;
    std::uint64_t random_state = 0x6669656c6463686bULL;
    std::array<Block, blocks> source;
    for (Block& value : source) {
        value = random_block(random_state);
    }
    std::array<std::uint64_t, 2> coefficients;
    for (std::uint64_t& coefficient : coefficients) {
        coefficient = splitmix64(random_state);
        coefficient |= static_cast<std::uint64_t>(coefficient == 0);
    }

    auto direct = source;
    auto bitsliced = source;
    auto rectangular = source;
    auto rectangular_avx2 = source;
    apply_checkpoint_chain<CheckpointImplementation::direct>(
        direct.data(), direct.size(), coefficients.data());
    apply_checkpoint_chain<CheckpointImplementation::bitsliced>(
        bitsliced.data(), bitsliced.size(), coefficients.data());
    apply_checkpoint_chain<CheckpointImplementation::rectangular>(
        rectangular.data(), rectangular.size(), coefficients.data());
    apply_checkpoint_chain<CheckpointImplementation::rectangular_avx2>(
        rectangular_avx2.data(), rectangular_avx2.size(), coefficients.data());
    for (std::size_t index = 0; index < blocks; ++index) {
        if (!equal_block(direct[index], bitsliced[index])) {
            throw std::runtime_error(
                "checkpoint direct and bitsliced evaluators disagree");
        }
        if (!equal_block(direct[index], rectangular[index])) {
            throw std::runtime_error(
                "checkpoint direct and rectangular evaluators disagree");
        }
        if (!equal_block(direct[index], rectangular_avx2[index])) {
            throw std::runtime_error(
                "checkpoint direct and rectangular AVX2 evaluators disagree");
        }
    }
}

void benchmark_checkpoint(
    CheckpointImplementation implementation,
    std::size_t word_blocks,
    unsigned trials,
    std::uint64_t seed) {
    if (word_blocks % checkpoint_epoch_blocks != 0) {
        throw std::invalid_argument("word-blocks must be divisible by 256");
    }
    const std::size_t epochs = word_blocks / checkpoint_epoch_blocks;
    std::uint64_t random_state = seed;
    std::vector<Block> source(word_blocks);
    std::vector<Block> working(word_blocks);
    for (Block& value : source) {
        value = random_block(random_state);
    }
    std::vector<std::uint64_t> coefficients(epochs - 1);
    for (std::uint64_t& coefficient : coefficients) {
        coefficient = splitmix64(random_state);
        coefficient |= static_cast<std::uint64_t>(coefficient == 0);
    }

    std::memcpy(working.data(), source.data(), source.size() * sizeof(Block));
    invoke_checkpoint(
        implementation, working.data(), word_blocks, coefficients.data());

    std::vector<double> milliseconds;
    milliseconds.reserve(trials);
    std::uint64_t final_checksum = 0;
    for (unsigned trial = 0; trial < trials; ++trial) {
        std::memcpy(working.data(), source.data(), source.size() * sizeof(Block));
        const auto begin = std::chrono::steady_clock::now();
        invoke_checkpoint(
            implementation, working.data(), word_blocks, coefficients.data());
        const auto end = std::chrono::steady_clock::now();
        milliseconds.push_back(
            std::chrono::duration<double, std::milli>(end - begin).count());
        final_checksum ^= std::rotl(checksum(working), trial & 63);
    }
    std::sort(milliseconds.begin(), milliseconds.end());
    const double median = milliseconds[milliseconds.size() / 2];
    const double minimum = milliseconds.front();
    const double gib = static_cast<double>(word_blocks * sizeof(Block)) /
        static_cast<double>(std::uint64_t{1} << 30);

    std::cout << std::fixed << std::setprecision(6)
              << "{\n"
              << "  \"schema\": \"riffle-transposed-field-checkpoint-benchmark-v1\",\n"
              << "  \"implementation\": \""
              << implementation_name(implementation) << "\",\n"
              << "  \"state_blocks\": " << checkpoint_state_blocks << ",\n"
              << "  \"epoch_blocks\": " << checkpoint_epoch_blocks << ",\n"
              << "  \"word_blocks\": " << word_blocks << ",\n"
              << "  \"epochs\": " << epochs << ",\n"
              << "  \"field_maps\": " << epochs - 1 << ",\n"
              << "  \"trials\": " << trials << ",\n"
              << "  \"median_ms\": " << median << ",\n"
              << "  \"minimum_ms\": " << minimum << ",\n"
              << "  \"median_input_gib_per_second\": "
              << gib / (median / 1000.0) << ",\n"
              << "  \"checksum\": \"0x" << std::hex << final_checksum
              << std::dec << "\",\n"
              << "  \"scope\": \"Complete reverse-order checkpoint inner on 128-bit elements; setup and input reset are outside timing.\"\n"
              << "}\n";
}

void benchmark_maps(
    CheckpointImplementation implementation,
    std::size_t map_count,
    unsigned trials,
    std::uint64_t seed) {
    if (implementation != CheckpointImplementation::bitsliced &&
        implementation != CheckpointImplementation::rectangular_avx2) {
        throw std::invalid_argument(
            "map-only benchmark supports bitsliced or rectangular-avx2");
    }
    std::uint64_t random_state = seed;
    std::vector<std::uint64_t> coefficients(map_count);
    for (std::uint64_t& coefficient : coefficients) {
        coefficient = splitmix64(random_state);
        coefficient |= static_cast<std::uint64_t>(coefficient == 0);
    }
    alignas(64) std::array<Block, 128> initial{};
    for (unsigned row = 0; row < checkpoint_state_blocks; ++row) {
        initial[row] = random_block(random_state);
    }
    alignas(64) std::array<Block, 128> working;

    auto invoke = [&] {
        if (implementation == CheckpointImplementation::bitsliced) {
            for (std::size_t map = 0; map < map_count; ++map) {
                transform_field_square<checkpoint_state_blocks>(
                    working.data(), coefficients[map]);
            }
        } else {
            for (std::size_t map = 0; map < map_count; ++map) {
                transform_field_rectangle64<transpose64_standard_avx2>(
                    working.data(), coefficients[map]);
            }
        }
    };

    working = initial;
    invoke();
    std::vector<double> milliseconds;
    milliseconds.reserve(trials);
    std::uint64_t final_checksum = 0;
    for (unsigned trial = 0; trial < trials; ++trial) {
        working = initial;
        const auto begin = std::chrono::steady_clock::now();
        invoke();
        const auto end = std::chrono::steady_clock::now();
        milliseconds.push_back(
            std::chrono::duration<double, std::milli>(end - begin).count());
        std::uint64_t trial_checksum = 0;
        for (unsigned row = 0; row < checkpoint_state_blocks; ++row) {
            trial_checksum ^= static_cast<std::uint64_t>(
                _mm_cvtsi128_si64(working[row].value));
            trial_checksum = std::rotl(trial_checksum, 7);
            trial_checksum ^= static_cast<std::uint64_t>(
                _mm_extract_epi64(working[row].value, 1));
        }
        final_checksum ^= std::rotl(trial_checksum, trial & 63);
    }
    std::sort(milliseconds.begin(), milliseconds.end());
    const double median = milliseconds[milliseconds.size() / 2];
    const double minimum = milliseconds.front();
    std::cout << std::fixed << std::setprecision(6)
              << "{\n"
              << "  \"schema\": \"riffle-field-checkpoint-map-only-benchmark-v1\",\n"
              << "  \"implementation\": \""
              << implementation_name(implementation) << "\",\n"
              << "  \"map_count\": " << map_count << ",\n"
              << "  \"trials\": " << trials << ",\n"
              << "  \"median_ms\": " << median << ",\n"
              << "  \"minimum_ms\": " << minimum << ",\n"
              << "  \"median_nanoseconds_per_map\": "
              << median * 1.0e6 / static_cast<double>(map_count) << ",\n"
              << "  \"minimum_nanoseconds_per_map\": "
              << minimum * 1.0e6 / static_cast<double>(map_count) << ",\n"
              << "  \"checksum\": \"0x" << std::hex << final_checksum
              << std::dec << "\"\n"
              << "}\n";
}

} // namespace

int main(int argc, char** argv) {
    try {
        CheckpointImplementation implementation =
            CheckpointImplementation::bitsliced;
        std::size_t word_blocks = std::size_t{1} << 21;
        unsigned trials = 15;
        std::uint64_t seed = 0x726966666c652d63ULL;
        bool self_test = false;
        bool map_only = false;
        std::size_t map_count = std::size_t{1} << 18;
        for (int index = 1; index < argc; ++index) {
            const std::string option = argv[index];
            auto value = [&]() -> std::string {
                if (++index >= argc) {
                    throw std::invalid_argument("missing value for " + option);
                }
                return argv[index];
            };
            if (option == "--implementation") {
                implementation = parse_checkpoint_implementation(value());
            } else if (option == "--word-blocks") {
                word_blocks = static_cast<std::size_t>(std::stoull(value()));
            } else if (option == "--trials") {
                trials = static_cast<unsigned>(std::stoul(value()));
            } else if (option == "--seed") {
                seed = std::stoull(value(), nullptr, 0);
            } else if (option == "--self-test") {
                self_test = true;
            } else if (option == "--map-only") {
                map_only = true;
            } else if (option == "--map-count") {
                map_count = static_cast<std::size_t>(std::stoull(value()));
            } else {
                throw std::invalid_argument("unknown option " + option);
            }
        }
        if (trials == 0) {
            throw std::invalid_argument("trials must be positive");
        }

        int cpu_features[4]{};
        __cpuidex(cpu_features, 7, 0);
        if ((self_test || implementation == CheckpointImplementation::bitsliced ||
             implementation == CheckpointImplementation::rectangular ||
             implementation == CheckpointImplementation::rectangular_avx2) &&
            (cpu_features[2] & (1 << 10)) == 0) {
            throw std::runtime_error("bitsliced kernel requires VPCLMULQDQ");
        }
        SetThreadAffinityMask(GetCurrentThread(), 1);
        SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_HIGHEST);
        if (self_test) {
            self_test_checkpoint();
            std::cout << "self_test,pass\n";
            return 0;
        }
        if (map_only) {
            benchmark_maps(implementation, map_count, trials, seed);
            return 0;
        }
        benchmark_checkpoint(implementation, word_blocks, trials, seed);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error," << error.what() << '\n';
        return 1;
    }
}
