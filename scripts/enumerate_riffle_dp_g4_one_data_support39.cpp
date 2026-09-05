#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{
    using u8 = std::uint8_t;
    using u16 = std::uint16_t;
    using u64 = std::uint64_t;

    constexpr u64 FieldReduction = 0x1b;
    constexpr u64 LocalGenerator = 0xF4845518B9582A1FULL;
    constexpr unsigned DataSymbols = 1 << 14;

    struct Hash
    {
        std::size_t operator()(u64 value) const
        {
            value ^= value >> 30;
            value *= 0xbf58476d1ce4e5b9ULL;
            value ^= value >> 27;
            value *= 0x94d049bb133111ebULL;
            return static_cast<std::size_t>(value ^ (value >> 31));
        }
    };

    u64 multiplyX(u64 value)
    {
        const auto high = value >> 63;
        return (value << 1) ^ (FieldReduction & (u64{0} - high));
    }

    unsigned localSupport(u64 message)
    {
        u64 low = 0;
        u64 high = 0;
        while (message)
        {
            const auto bit = std::countr_zero(message);
            low ^= LocalGenerator << bit;
            high ^=
                (bit == 0 ? 0 : LocalGenerator >> (64 - bit)) |
                (u64{1} << 63);
            message &= message - 1;
        }
        unsigned support = 0;
        for (unsigned slot = 0; slot < 16; ++slot)
        {
            support += ((low >> (4 * slot)) & 0xf) != 0;
            support += ((high >> (4 * slot)) & 0xf) != 0;
        }
        return support;
    }

    struct Record
    {
        u16 exponent;
        u64 parameter;
        u8 totalSupport;
    };
}

int main(int argc, char** argv)
{
    try
    {
        if (argc != 4)
            throw std::invalid_argument(
                "usage: enumerate_riffle_dp_g4_one_data_support39 LOW.bin OUT.bin OUT.json");
        const std::string inputPath = argv[1];
        const std::string recordPath = argv[2];
        const std::string receiptPath = argv[3];

        std::ifstream input(inputPath, std::ios::binary);
        if (!input)
            throw std::runtime_error("one-data enumeration: could not open input");
        std::unordered_map<u64, u8, Hash> supports;
        supports.reserve(80'000);
        std::vector<std::pair<u64, u8>> parameters;
        while (true)
        {
            u8 support = 0;
            u64 message = 0;
            input.read(reinterpret_cast<char*>(&support), 1);
            if (!input)
            {
                if (input.eof())
                    break;
                throw std::runtime_error("one-data enumeration: support read failed");
            }
            input.read(reinterpret_cast<char*>(&message), sizeof(message));
            if (!input || support < 11 || support > 13)
                throw std::runtime_error("one-data enumeration: invalid input record");
            if (localSupport(message) != support)
                throw std::runtime_error("one-data enumeration: local re-encoding failed");
            if (!supports.emplace(message, support).second)
                throw std::runtime_error("one-data enumeration: duplicate input message");
            parameters.emplace_back(message, support);
        }
        if (parameters.size() != 38'560)
            throw std::runtime_error("one-data enumeration: incomplete support-13 input");

        std::ofstream records(recordPath, std::ios::binary);
        if (!records)
            throw std::runtime_error("one-data enumeration: could not open record output");
        std::array<u64, 40> counts{};
        u64 tested = 0;
        for (const auto& [parameter, parameterSupport] : parameters)
        {
            u64 weighted = parameter;
            for (unsigned exponent = 0; exponent < DataSymbols; ++exponent)
            {
                ++tested;
                const auto found = supports.find(weighted);
                if (found != supports.end())
                {
                    const auto totalSupport = static_cast<unsigned>(2 * parameterSupport + found->second);
                    if (totalSupport > 39)
                        throw std::runtime_error("one-data enumeration: support bound failed");
                    ++counts[totalSupport];
                    const Record record{
                        static_cast<u16>(exponent),
                        parameter,
                        static_cast<u8>(totalSupport)};
                    records.write(reinterpret_cast<const char*>(&record.exponent), sizeof(record.exponent));
                    records.write(reinterpret_cast<const char*>(&record.parameter), sizeof(record.parameter));
                    records.write(reinterpret_cast<const char*>(&record.totalSupport), sizeof(record.totalSupport));
                }
                weighted = multiplyX(weighted);
            }
        }
        if (!records)
            throw std::runtime_error("one-data enumeration: record write failed");

        std::ofstream receipt(receiptPath);
        receipt << "{\n"
            << "  \"schema\": \"riffle-dp-g4-one-data-support39-enumeration-v1\",\n"
            << "  \"candidate\": \"Riffle DP g=4@g0-v1\",\n"
            << "  \"evidence_label\": \"EXACT\",\n"
            << "  \"local_message_count\": " << parameters.size() << ",\n"
            << "  \"coefficient_exponents\": " << DataSymbols << ",\n"
            << "  \"tested_pairs\": " << tested << ",\n"
            << "  \"outer_word_counts\": {\n";
        bool first = true;
        u64 total = 0;
        for (unsigned support = 0; support < counts.size(); ++support)
        {
            if (!counts[support])
                continue;
            if (!first)
                receipt << ",\n";
            first = false;
            receipt << "    \"" << support << "\": " << counts[support];
            total += counts[support];
        }
        receipt << "\n  },\n"
            << "  \"record_count\": " << total << ",\n"
            << "  \"record_format\": \"little-endian uint16 exponent, uint64 parameter, uint8 total support\",\n"
            << "  \"scope\": \"All one-data words (t,t,x^i t) whose three local words have support at most 13.\"\n"
            << "}\n";
        if (!receipt)
            throw std::runtime_error("one-data enumeration: receipt write failed");

        std::cout << "candidate=Riffle DP g=4@g0-v1\n";
        std::cout << "tested_pairs=" << tested << '\n';
        for (unsigned support = 0; support < counts.size(); ++support)
            if (counts[support])
                std::cout << "support_" << support << "=" << counts[support] << '\n';
        std::cout << "records=" << total << '\n';
        std::cout << "status=EXACT_ONE_DATA_SUPPORT39_ENUMERATION\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
