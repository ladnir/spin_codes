#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#if defined(_MSC_VER)
#include <intrin.h>
#endif

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace
{
    constexpr std::size_t outerBits = 250;
    constexpr std::size_t outerWords = 4;
    constexpr std::size_t outerRows = 8576;
    constexpr std::size_t fanoutLayers = 56;
    constexpr std::size_t sourceBits = 31;
    constexpr std::size_t targetBits = 33;
    constexpr std::uint64_t finalWordMask = (std::uint64_t{1} << 58) - 1;

    struct PackedRow
    {
        std::array<std::uint64_t, outerWords> words{};
    };

    static_assert(sizeof(PackedRow) == 32);

    struct alignas(64) FanoutLayer
    {
        std::array<std::uint64_t, outerWords> sources{};
        std::array<std::uint64_t, outerWords> targets{};
    };

    static_assert(sizeof(FanoutLayer) == 64);

    volatile std::uint64_t benchmarkSink = 0;

    std::string cpuBrand()
    {
#if defined(_MSC_VER) && (defined(_M_IX86) || defined(_M_X64))
        std::array<int, 4> registers{};
        __cpuid(registers.data(), 0x80000000);
        if (static_cast<unsigned>(registers[0]) < 0x80000004u)
            return "unavailable";

        std::array<char, 49> brand{};
        for (int leaf = 0; leaf < 3; ++leaf)
        {
            __cpuid(registers.data(), 0x80000002 + leaf);
            std::memcpy(brand.data() + 16 * leaf, registers.data(), 16);
        }
        const std::string value(brand.data());
        const auto first = value.find_first_not_of(' ');
        const auto last = value.find_last_not_of(' ');
        return first == std::string::npos
            ? "unavailable"
            : value.substr(first, last - first + 1);
#else
        return "unavailable";
#endif
    }

    std::string compilerName()
    {
#if defined(_MSC_VER)
        return "MSVC " + std::to_string(_MSC_VER);
#elif defined(__clang__)
        return "Clang " __clang_version__;
#elif defined(__GNUC__)
        return "GCC " __VERSION__;
#else
        return "unavailable";
#endif
    }

    std::string jsonEscape(const std::string& value)
    {
        std::string escaped;
        escaped.reserve(value.size());
        for (const char character : value)
        {
            if (character == '\\' || character == '"')
                escaped.push_back('\\');
            escaped.push_back(character);
        }
        return escaped;
    }

    std::uint64_t splitmix64(std::uint64_t& state) noexcept
    {
        state += 0x9e3779b97f4a7c15ULL;
        std::uint64_t value = state;
        value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
        value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
        return value ^ (value >> 31);
    }

    std::uint64_t boundedRandom(std::uint64_t& state, std::uint64_t bound)
    {
        if (bound == 0)
            throw std::invalid_argument("zero random bound");
        const std::uint64_t threshold = (std::uint64_t{0} - bound) % bound;
        for (;;)
        {
            const std::uint64_t value = splitmix64(state);
            if (value >= threshold)
                return value % bound;
        }
    }

    void setBit(std::array<std::uint64_t, outerWords>& mask, std::size_t bit)
    {
        mask[bit >> 6] |= std::uint64_t{1} << (bit & 63);
    }

    std::vector<FanoutLayer> buildSchedules(std::uint64_t seed)
    {
        std::vector<FanoutLayer> schedules(outerRows * fanoutLayers);
        std::array<std::uint8_t, outerBits> coordinates{};
        for (std::size_t row = 0; row < outerRows; ++row)
        {
            for (std::size_t layer = 0; layer < fanoutLayers; ++layer)
            {
                std::iota(coordinates.begin(), coordinates.end(), std::uint8_t{0});
                for (std::size_t pick = 0; pick < sourceBits + targetBits; ++pick)
                {
                    const std::size_t remaining = outerBits - pick;
                    const std::size_t offset = static_cast<std::size_t>(
                        boundedRandom(seed, remaining));
                    std::swap(coordinates[pick], coordinates[pick + offset]);
                }

                auto& schedule = schedules[row * fanoutLayers + layer];
                for (std::size_t index = 0; index < sourceBits; ++index)
                    setBit(schedule.sources, coordinates[index]);
                for (std::size_t index = 0; index < targetBits; ++index)
                    setBit(schedule.targets, coordinates[sourceBits + index]);
            }
        }
        return schedules;
    }

    void validateSchedules(const std::vector<FanoutLayer>& schedules)
    {
        if (schedules.size() != outerRows * fanoutLayers)
            throw std::runtime_error("schedule count mismatch");
        for (const auto& schedule : schedules)
        {
            unsigned sources = 0;
            unsigned targets = 0;
            for (std::size_t word = 0; word < outerWords; ++word)
            {
                if (schedule.sources[word] & schedule.targets[word])
                    throw std::runtime_error("source and target masks overlap");
                sources += std::popcount(schedule.sources[word]);
                targets += std::popcount(schedule.targets[word]);
            }
            if (sources != sourceBits || targets != targetBits)
                throw std::runtime_error("schedule weight mismatch");
            if ((schedule.sources.back() | schedule.targets.back()) & ~finalWordMask)
                throw std::runtime_error("schedule uses a padding bit");
        }
    }

    inline void applyLayer(PackedRow& row, const FanoutLayer& layer) noexcept
    {
        const std::uint64_t folded =
            (row.words[0] & layer.sources[0]) ^
            (row.words[1] & layer.sources[1]) ^
            (row.words[2] & layer.sources[2]) ^
            (row.words[3] & layer.sources[3]);
        const std::uint64_t mask = std::uint64_t{0} - (std::popcount(folded) & 1u);
        row.words[0] ^= layer.targets[0] & mask;
        row.words[1] ^= layer.targets[1] & mask;
        row.words[2] ^= layer.targets[2] & mask;
        row.words[3] ^= layer.targets[3] & mask;
    }

    void applyFanout56(
        std::vector<PackedRow>& rows,
        const std::vector<FanoutLayer>& schedules) noexcept
    {
        for (std::size_t row = 0; row < outerRows; ++row)
        {
            PackedRow& value = rows[row];
            const FanoutLayer* layers = schedules.data() + row * fanoutLayers;
            for (std::size_t layer = 0; layer < fanoutLayers; ++layer)
                applyLayer(value, layers[layer]);
        }
    }

    void invertFanout56(
        std::vector<PackedRow>& rows,
        const std::vector<FanoutLayer>& schedules) noexcept
    {
        for (std::size_t row = 0; row < outerRows; ++row)
        {
            PackedRow& value = rows[row];
            const FanoutLayer* layers = schedules.data() + row * fanoutLayers;
            for (std::size_t layer = fanoutLayers; layer-- > 0;)
                applyLayer(value, layers[layer]);
        }
    }

    bool bit(const PackedRow& row, std::size_t index) noexcept
    {
        return ((row.words[index >> 6] >> (index & 63)) & 1) != 0;
    }

    void toggle(PackedRow& row, std::size_t index) noexcept
    {
        row.words[index >> 6] ^= std::uint64_t{1} << (index & 63);
    }

    void applyLayerScalar(PackedRow& row, const FanoutLayer& layer)
    {
        bool parity = false;
        for (std::size_t index = 0; index < outerBits; ++index)
            if (bit(PackedRow{layer.sources}, index))
                parity ^= bit(row, index);
        if (!parity)
            return;
        for (std::size_t index = 0; index < outerBits; ++index)
            if (bit(PackedRow{layer.targets}, index))
                toggle(row, index);
    }

    void selfTest(
        const std::vector<PackedRow>& input,
        const std::vector<FanoutLayer>& schedules)
    {
        for (std::size_t row = 0; row < 32; ++row)
        {
            PackedRow packed = input[row];
            PackedRow scalar = input[row];
            const FanoutLayer* layers = schedules.data() + row * fanoutLayers;
            for (std::size_t layer = 0; layer < fanoutLayers; ++layer)
            {
                applyLayer(packed, layers[layer]);
                applyLayerScalar(scalar, layers[layer]);
            }
            if (packed.words != scalar.words)
                throw std::runtime_error("packed and scalar fanout disagree");
        }

        std::vector<PackedRow> roundTrip = input;
        applyFanout56(roundTrip, schedules);
        invertFanout56(roundTrip, schedules);
        if (std::memcmp(roundTrip.data(), input.data(), input.size() * sizeof(PackedRow)))
            throw std::runtime_error("reverse fanout did not recover the input");
    }

    std::uint64_t checksum(const std::vector<PackedRow>& rows) noexcept
    {
        std::uint64_t value = 0x243f6a8885a308d3ULL;
        for (const auto& row : rows)
            for (const auto word : row.words)
                value = std::rotl(value, 9) ^ word;
        return value;
    }

    double milliseconds(const std::chrono::steady_clock::duration duration)
    {
        return std::chrono::duration<double, std::milli>(duration).count();
    }

    template<typename Operation>
    double timeBatch(Operation&& operation, std::size_t batch)
    {
        const auto start = std::chrono::steady_clock::now();
        for (std::size_t iteration = 0; iteration < batch; ++iteration)
            operation();
        const auto stop = std::chrono::steady_clock::now();
        return milliseconds(stop - start) / static_cast<double>(batch);
    }

    double quantile(std::vector<double> values, double fraction)
    {
        std::sort(values.begin(), values.end());
        const double position = fraction * static_cast<double>(values.size() - 1);
        const std::size_t lower = static_cast<std::size_t>(position);
        const std::size_t upper = std::min(lower + 1, values.size() - 1);
        const double share = position - static_cast<double>(lower);
        return values[lower] * (1.0 - share) + values[upper] * share;
    }

    void pinThread(unsigned cpu)
    {
#if defined(_WIN32)
        if (cpu >= 64)
            throw std::invalid_argument("CPU index exceeds one affinity group");
        const DWORD_PTR mask = DWORD_PTR{1} << cpu;
        if (!SetThreadAffinityMask(GetCurrentThread(), mask))
            throw std::runtime_error("SetThreadAffinityMask failed");
#else
        (void)cpu;
#endif
    }

    struct Options
    {
        std::size_t trials = 31;
        std::size_t warmup = 4;
        std::size_t batch = 32;
        unsigned cpu = 2;
        std::string output;
    };

    Options parseOptions(int argc, char** argv)
    {
        Options result;
        for (int index = 1; index < argc; ++index)
        {
            const std::string argument = argv[index];
            auto value = [&]() -> std::string {
                if (++index >= argc)
                    throw std::invalid_argument("missing option value");
                return argv[index];
            };
            if (argument == "--trials")
                result.trials = std::stoull(value());
            else if (argument == "--warmup")
                result.warmup = std::stoull(value());
            else if (argument == "--batch")
                result.batch = std::stoull(value());
            else if (argument == "--cpu")
                result.cpu = static_cast<unsigned>(std::stoul(value()));
            else if (argument == "--output")
                result.output = value();
            else
                throw std::invalid_argument("unknown option: " + argument);
        }
        if (result.trials < 3 || result.batch == 0 || result.output.empty())
            throw std::invalid_argument("invalid benchmark options");
        return result;
    }

    void printSamples(std::ostream& output, const std::vector<double>& values)
    {
        output << '[';
        for (std::size_t index = 0; index < values.size(); ++index)
        {
            if (index)
                output << ',';
            output << std::setprecision(12) << values[index];
        }
        output << ']';
    }
}

int main(int argc, char** argv)
{
    try
    {
        const Options options = parseOptions(argc, argv);
        pinThread(options.cpu);

        std::uint64_t inputSeed = 0xb250124056fa901dULL;
        std::vector<PackedRow> input(outerRows);
        for (auto& row : input)
        {
            for (auto& word : row.words)
                word = splitmix64(inputSeed);
            row.words.back() &= finalWordMask;
        }
        std::vector<PackedRow> working(outerRows);

        const auto setupStart = std::chrono::steady_clock::now();
        const auto schedules = buildSchedules(0x31'33'00'56'b250'0124ULL);
        const auto setupStop = std::chrono::steady_clock::now();
        const double setupMs = milliseconds(setupStop - setupStart);
        validateSchedules(schedules);
        selfTest(input, schedules);

        auto zero = [&] {
            std::memcpy(working.data(), input.data(), input.size() * sizeof(PackedRow));
        };
        auto layer56 = [&] {
            std::memcpy(working.data(), input.data(), input.size() * sizeof(PackedRow));
            applyFanout56(working, schedules);
        };

        for (std::size_t iteration = 0; iteration < options.warmup; ++iteration)
        {
            zero();
            benchmarkSink ^= checksum(working);
            layer56();
            benchmarkSink ^= checksum(working);
        }

        std::vector<double> zeroSamples;
        std::vector<double> layer56Samples;
        zeroSamples.reserve(options.trials);
        layer56Samples.reserve(options.trials);
        for (std::size_t trial = 0; trial < options.trials; ++trial)
        {
            if ((trial & 1) == 0)
            {
                zeroSamples.push_back(timeBatch(zero, options.batch));
                benchmarkSink ^= checksum(working);
                layer56Samples.push_back(timeBatch(layer56, options.batch));
                benchmarkSink ^= checksum(working);
            }
            else
            {
                layer56Samples.push_back(timeBatch(layer56, options.batch));
                benchmarkSink ^= checksum(working);
                zeroSamples.push_back(timeBatch(zero, options.batch));
                benchmarkSink ^= checksum(working);
            }
        }

        const double zeroMedian = quantile(zeroSamples, 0.5);
        const double layer56Median = quantile(layer56Samples, 0.5);
        const double incremental = layer56Median - zeroMedian;
        const double logicalRowsMiB =
            static_cast<double>(outerRows * outerBits) / 8.0 / (1024.0 * 1024.0);
        const double packedRowsMiB =
            static_cast<double>(input.size() * sizeof(PackedRow)) /
            (1024.0 * 1024.0);
        const std::size_t scheduleBytes = schedules.size() * sizeof(FanoutLayer);
        const double scheduleMiB =
            static_cast<double>(scheduleBytes) /
            (1024.0 * 1024.0);
        const double scheduleStreamGbps =
            static_cast<double>(scheduleBytes) / incremental / 1'000'000.0;
        const std::uint64_t inputChecksum = checksum(input);
        zero();
        const std::uint64_t zeroChecksum = checksum(working);
        layer56();
        const std::uint64_t layer56Checksum = checksum(working);

        std::ofstream output(options.output, std::ios::binary);
        if (!output)
            throw std::runtime_error("failed to open benchmark output");
        output << "{\n";
        output << "  \"schema\": \"bch250-rowlocal-fanout56-benchmark-v1\",\n";
        output << "  \"status\": \"MEASURED_WRAPPER_ONLY\",\n";
        output << "  \"environment\": {\n";
        output << "    \"cpu_brand\": \"" << jsonEscape(cpuBrand()) << "\",\n";
        output << "    \"compiler\": \"" << jsonEscape(compilerName()) << "\",\n";
        output << "    \"pointer_bits\": " << 8 * sizeof(void*) << "\n";
        output << "  },\n";
        output << "  \"parameters\": {\n";
        output << "    \"outer_bits\": " << outerBits << ",\n";
        output << "    \"outer_rows\": " << outerRows << ",\n";
        output << "    \"fanout_layers\": " << fanoutLayers << ",\n";
        output << "    \"source_bits_per_layer\": " << sourceBits << ",\n";
        output << "    \"target_bits_per_layer\": " << targetBits << ",\n";
        output << "    \"trials\": " << options.trials << ",\n";
        output << "    \"warmup_pairs\": " << options.warmup << ",\n";
        output << "    \"batch_encodings_per_sample\": " << options.batch << ",\n";
        output << "    \"pinned_logical_cpu\": " << options.cpu << "\n";
        output << "  },\n";
        output << "  \"correctness\": {\n";
        output << "    \"packed_matches_scalar_first_32_rows\": true,\n";
        output << "    \"reverse_layer_order_recovers_all_rows\": true,\n";
        output << "    \"source_target_disjoint_and_weights_valid\": true\n";
        output << "  },\n";
        output << "  \"memory\": {\n";
        output << "    \"logical_rows_mib\": " << std::setprecision(12)
               << logicalRowsMiB << ",\n";
        output << "    \"packed_rows_storage_mib\": " << packedRowsMiB << ",\n";
        output << "    \"schedule_mib\": " << scheduleMiB << ",\n";
        output << "    \"schedule_bytes\": " << scheduleBytes << ",\n";
        output << "    \"bytes_per_layer_schedule\": " << sizeof(FanoutLayer) << "\n";
        output << "  },\n";
        output << "  \"setup_ms\": " << setupMs << ",\n";
        output << "  \"zero_layers\": {\n";
        output << "    \"median_ms\": " << zeroMedian << ",\n";
        output << "    \"p10_ms\": " << quantile(zeroSamples, 0.1) << ",\n";
        output << "    \"p90_ms\": " << quantile(zeroSamples, 0.9) << ",\n";
        output << "    \"samples_ms\": ";
        printSamples(output, zeroSamples);
        output << "\n  },\n";
        output << "  \"fifty_six_layers\": {\n";
        output << "    \"median_ms\": " << layer56Median << ",\n";
        output << "    \"p10_ms\": " << quantile(layer56Samples, 0.1) << ",\n";
        output << "    \"p90_ms\": " << quantile(layer56Samples, 0.9) << ",\n";
        output << "    \"samples_ms\": ";
        printSamples(output, layer56Samples);
        output << "\n  },\n";
        output << "  \"comparison\": {\n";
        output << "    \"incremental_median_ms\": " << incremental << ",\n";
        output << "    \"ratio_56_over_0\": " << layer56Median / zeroMedian << ",\n";
        output << "    \"million_layer_applications_per_second\": "
               << static_cast<double>(outerRows * fanoutLayers) /
                      incremental / 1000.0 << ",\n";
        output << "    \"schedule_stream_gb_per_second\": "
               << scheduleStreamGbps << "\n";
        output << "  },\n";
        output << "  \"checksums_hex\": {\n";
        output << "    \"input\": \"0x" << std::hex << inputChecksum << "\",\n";
        output << "    \"zero_layers\": \"0x" << zeroChecksum << "\",\n";
        output << "    \"fifty_six_layers\": \"0x" << layer56Checksum << "\"\n";
        output << std::dec << "  },\n";
        output << "  \"scope\": [\n";
        output << "    \"Both timings copy the same 8576 packed 250-bit rows; the difference isolates the 56-layer mask kernel.\",\n";
        output << "    \"The benchmark starts from post-BCH rows and omits the local coordinate permutations, transpose, region permutations, and RM2Sub-S19.\",\n";
        output << "    \"Each row has an independent 56-layer schedule stored as source and target masks.\"\n";
        output << "  ]\n";
        output << "}\n";
        output.close();

        std::cout << std::fixed << std::setprecision(6)
                  << "setup_ms=" << setupMs << '\n'
                  << "zero_median_ms=" << zeroMedian << '\n'
                  << "fanout56_median_ms=" << layer56Median << '\n'
                  << "incremental_median_ms=" << incremental << '\n'
                  << "schedule_mib=" << scheduleMiB << '\n'
                  << "output=" << options.output << '\n';
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
