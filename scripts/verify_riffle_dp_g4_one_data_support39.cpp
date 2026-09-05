#include <bit>
#include <cstdint>
#include <fstream>
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
}

int main(int argc, char** argv)
{
    try
    {
        if (argc != 3)
            throw std::invalid_argument(
                "usage: verify_riffle_dp_g4_one_data_support39 LOW.bin RECORDS.bin");
        std::ifstream input(argv[1], std::ios::binary);
        std::ifstream records(argv[2], std::ios::binary);
        if (!input || !records)
            throw std::runtime_error("one-data verifier: input open failed");

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
                throw std::runtime_error("one-data verifier: local support read failed");
            }
            input.read(reinterpret_cast<char*>(&message), sizeof(message));
            if (!input || support < 11 || support > 13 || localSupport(message) != support)
                throw std::runtime_error("one-data verifier: invalid local message");
            if (!supports.emplace(message, support).second)
                throw std::runtime_error("one-data verifier: duplicate local message");
            parameters.emplace_back(message, support);
        }
        if (parameters.size() != 38'560)
            throw std::runtime_error("one-data verifier: local corpus size mismatch");

        u64 matched = 0;
        u64 tested = 0;
        for (const auto& [parameter, parameterSupport] : parameters)
        {
            auto weighted = parameter;
            for (unsigned exponent = 0; exponent < DataSymbols; ++exponent)
            {
                ++tested;
                const auto found = supports.find(weighted);
                if (found != supports.end())
                {
                    u16 recordExponent = 0;
                    u64 recordParameter = 0;
                    u8 recordSupport = 0;
                    records.read(reinterpret_cast<char*>(&recordExponent), sizeof(recordExponent));
                    records.read(reinterpret_cast<char*>(&recordParameter), sizeof(recordParameter));
                    records.read(reinterpret_cast<char*>(&recordSupport), sizeof(recordSupport));
                    const auto expectedSupport = static_cast<u8>(2 * parameterSupport + found->second);
                    if (!records || recordExponent != exponent ||
                        recordParameter != parameter || recordSupport != expectedSupport)
                        throw std::runtime_error("one-data verifier: record stream mismatch");
                    ++matched;
                }
                weighted = multiplyX(weighted);
            }
        }
        char extra = 0;
        records.read(&extra, 1);
        if (!records.eof())
            throw std::runtime_error("one-data verifier: trailing record data");
        if (tested != 631'767'040 || matched != 85'828)
            throw std::runtime_error("one-data verifier: frozen count mismatch");
        std::cout << "candidate=Riffle DP g=4@g0-v1\n";
        std::cout << "tested_pairs=" << tested << '\n';
        std::cout << "matched_records=" << matched << '\n';
        std::cout << "status=EXACT_ONE_DATA_SUPPORT39_ENUMERATION_VERIFIED\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
