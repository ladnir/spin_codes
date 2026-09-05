#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace
{
	using u8 = std::uint8_t;
	using u16 = std::uint16_t;
	using u32 = std::uint32_t;
	using u64 = std::uint64_t;

	constexpr u64 dataSymbols = u64{ 1 } << 14;
	constexpr u64 fieldReduction = 0x1b;

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

	struct Candidate
	{
		std::array<u32, 3> words{};
	};

	struct AdditiveTriple
	{
		u64 left = 0;
		u64 right = 0;
		u64 third = 0;
		u8 totalSupport = 0;
	};

	u64 multiplyX(u64 value)
	{
		const auto high = value >> 63;
		return (value << 1) ^ (fieldReduction & (u64{ 0 } - high));
	}

	u64 multiply(u64 left, u64 right)
	{
		u64 result = 0;
		for (u64 bit = 0; bit < 64; ++bit)
		{
			result ^= left & (u64{ 0 } - (right & 1));
			left = multiplyX(left);
			right >>= 1;
		}
		return result;
	}

	u64 power(u64 value, u64 exponent)
	{
		u64 result = 1;
		while (exponent)
		{
			if (exponent & 1)
				result = multiply(result, value);
			value = multiply(value, value);
			exponent >>= 1;
		}
		return result;
	}

	u64 inverse(u64 value)
	{
		if (!value)
			throw std::invalid_argument("zero has no field inverse");
		return power(value, std::numeric_limits<u64>::max() - 1);
	}

	std::unordered_map<u64, u8, Hash> readMessages(
		const std::string& path,
		std::array<std::vector<u64>, 14>& bySupport)
	{
		std::ifstream input(path, std::ios::binary);
		if (!input)
			throw std::runtime_error("could not open low-message file");
		std::unordered_map<u64, u8, Hash> support;
		support.reserve(50'000);
		while (true)
		{
			u8 packetSupport = 0;
			u64 message = 0;
			input.read(reinterpret_cast<char*>(&packetSupport), 1);
			if (!input)
			{
				if (input.eof())
					break;
				throw std::runtime_error("failed reading support byte");
			}
			input.read(reinterpret_cast<char*>(&message), sizeof(message));
			if (!input || packetSupport >= bySupport.size())
				throw std::runtime_error("invalid low-message record");
			if (!support.emplace(message, packetSupport).second)
				throw std::runtime_error("duplicate low-message record");
			bySupport[packetSupport].push_back(message);
		}
		return support;
	}

	template<typename Function>
	void forEachAllowedOrderedPair(
		const std::array<std::vector<u64>, 14>& bySupport,
		Function&& function)
	{
		for (u64 leftSupport = 11; leftSupport <= 13; ++leftSupport)
		{
			for (u64 rightSupport = 11; rightSupport <= 13; ++rightSupport)
			{
				if (leftSupport + rightSupport + 11 > 35)
					continue;
				for (const auto left : bySupport[leftSupport])
				{
					for (const auto right : bySupport[rightSupport])
						function(left, static_cast<u8>(leftSupport), right,
							static_cast<u8>(rightSupport));
				}
			}
		}
	}
}

int main(int argc, char** argv)
{
	try
	{
		if (argc != 2)
			throw std::invalid_argument("usage: classify_riffle_dp_g4_low_triples FILE");
		std::array<std::vector<u64>, 14> bySupport;
		const auto support = readMessages(argv[1], bySupport);
		if (support.size() != 38'560 || bySupport[11].size() != 20 ||
			bySupport[12].size() != 1'526 || bySupport[13].size() != 37'014)
			throw std::runtime_error("low-message file is not the complete support-13 receipt");

		std::vector<u64> powers(dataSymbols);
		std::unordered_map<u64, u16, Hash> powerLog;
		powerLog.reserve(dataSymbols * 2);
		powers[0] = 1;
		for (u64 exponent = 1; exponent < dataSymbols; ++exponent)
			powers[exponent] = multiplyX(powers[exponent - 1]);
		for (u64 exponent = 0; exponent < dataSymbols; ++exponent)
			powerLog.emplace(powers[exponent], static_cast<u16>(exponent));
		if (powerLog.size() != dataSymbols)
			throw std::runtime_error("coefficient schedule is not distinct");

		std::unordered_map<u64, Candidate, Hash> repeatedDataCandidates;
		repeatedDataCandidates.reserve(1'000'000);
		for (u64 parameterSupport = 11; parameterSupport <= 12; ++parameterSupport)
		{
			const auto maximumThirdSupport = 35 - 2 * parameterSupport;
			for (const auto parameter : bySupport[parameterSupport])
			{
				const auto parameterInverse = inverse(parameter);
				for (u64 thirdSupport = 11; thirdSupport <= maximumThirdSupport; ++thirdSupport)
				{
					for (const auto third : bySupport[thirdSupport])
					{
						const auto coefficient = multiply(third, parameterInverse);
						auto& candidate = repeatedDataCandidates[coefficient];
						++candidate.words[2 * parameterSupport + thirdSupport - 33];
					}
				}
			}
		}

		std::array<u64, 3> oneWords{};
		std::array<u64, 3> oneDirty{};
		for (u64 index = 0; index < dataSymbols; ++index)
		{
			const auto found = repeatedDataCandidates.find(powers[index]);
			if (found == repeatedDataCandidates.end())
				continue;
			bool seen = false;
			for (u64 layer = 0; layer < 3; ++layer)
			{
				oneWords[layer] += found->second.words[layer];
				if (!seen && found->second.words[layer])
				{
					++oneDirty[layer];
					seen = true;
				}
			}
		}

		std::array<u64, 3> pairP1Words{};
		std::array<u64, 3> pairP1Dirty{};
		for (u64 left = 0; left + 1 < dataSymbols; ++left)
		{
			for (u64 right = left + 1; right < dataSymbols; ++right)
			{
				const auto found = repeatedDataCandidates.find(powers[left] ^ powers[right]);
				if (found == repeatedDataCandidates.end())
					continue;
				bool seen = false;
				for (u64 layer = 0; layer < 3; ++layer)
				{
					pairP1Words[layer] += found->second.words[layer];
					if (!seen && found->second.words[layer])
					{
						++pairP1Dirty[layer];
						seen = true;
					}
				}
			}
		}

		std::unordered_map<u64, u64, Hash> inverseCache;
		inverseCache.reserve(support.size());
		auto cachedInverse = [&](u64 value) {
			const auto found = inverseCache.find(value);
			if (found != inverseCache.end())
				return found->second;
			const auto result = inverse(value);
			inverseCache.emplace(value, result);
			return result;
		};

		std::vector<std::array<u64, 3>> pairP0ByGap(dataSymbols);
		std::vector<AdditiveTriple> additiveTriples;
		forEachAllowedOrderedPair(bySupport,
			[&](u64 left, u8 leftSupport, u64 right, u8 rightSupport) {
				const auto third = left ^ right;
				const auto foundThird = support.find(third);
				if (foundThird == support.end())
					return;
				const auto total = u64{ leftSupport } + rightSupport + foundThird->second;
				if (total > 35)
					return;

				const auto ratio = multiply(left, cachedInverse(right));
				const auto foundGap = powerLog.find(ratio);
				if (foundGap != powerLog.end() && foundGap->second != 0)
					++pairP0ByGap[foundGap->second][total - 33];

				additiveTriples.push_back(
					{ left, right, third, static_cast<u8>(total) });
			});

		std::array<u64, 3> pairP0Words{};
		std::array<u64, 3> pairP0Dirty{};
		for (u64 gap = 1; gap < dataSymbols; ++gap)
		{
			bool seen = false;
			for (u64 layer = 0; layer < 3; ++layer)
			{
				pairP0Words[layer] += pairP0ByGap[gap][layer] * (dataSymbols - gap);
				if (!seen && pairP0ByGap[gap][layer])
				{
					pairP0Dirty[layer] += dataSymbols - gap;
					seen = true;
				}
			}
		}

		std::unordered_map<u32, std::array<u64, 3>> tripleByGap;
		tripleByGap.reserve(1024);
		for (const auto& triple : additiveTriples)
		{
			std::unordered_map<u64, u16, Hash> scaledThirdLog;
			scaledThirdLog.reserve(dataSymbols * 2);
			auto scaledThird = triple.third;
			for (u64 gap = 1; gap < dataSymbols; ++gap)
			{
				scaledThird = multiplyX(scaledThird);
				scaledThirdLog.emplace(scaledThird, static_cast<u16>(gap));
			}

			auto scaledRight = triple.right;
			for (u64 gap1 = 1; gap1 + 1 < dataSymbols; ++gap1)
			{
				scaledRight = multiplyX(scaledRight);
				const auto target = triple.left ^ scaledRight;
				const auto foundGap2 = scaledThirdLog.find(target);
				if (foundGap2 == scaledThirdLog.end() || foundGap2->second <= gap1)
					continue;
				const auto key = static_cast<u32>(gap1 << 14) | foundGap2->second;
				++tripleByGap[key][triple.totalSupport - 33];
			}
		}

		std::array<u64, 3> tripleWords{};
		std::array<u64, 3> tripleDirty{};
		for (const auto& [key, counts] : tripleByGap)
		{
			const auto gap2 = key & ((u32{ 1 } << 14) - 1);
			bool seen = false;
			for (u64 layer = 0; layer < 3; ++layer)
			{
				tripleWords[layer] += counts[layer] * (dataSymbols - gap2);
				if (!seen && counts[layer])
				{
					tripleDirty[layer] += dataSymbols - gap2;
					seen = true;
				}
			}
		}

		std::cout << "candidate=Riffle DP g=4\n";
		std::cout << "local_messages=38560\n";
		for (u64 layer = 0; layer < 3; ++layer)
		{
			std::cout << "support=" << 33 + layer
				<< " one_words=" << oneWords[layer]
				<< " one_new_dirty_classes=" << oneDirty[layer]
				<< " pair_p0_words=" << pairP0Words[layer]
				<< " pair_p0_new_dirty_classes=" << pairP0Dirty[layer]
				<< " pair_p1_words=" << pairP1Words[layer]
				<< " pair_p1_new_dirty_classes=" << pairP1Dirty[layer]
				<< " triple_words=" << tripleWords[layer]
				<< " triple_new_dirty_classes=" << tripleDirty[layer]
				<< '\n';
		}
		std::cout << "triple_additive_ordered_candidates=" << additiveTriples.size() << '\n';
		std::cout << "triple_dirty_gap_classes=" << tripleByGap.size() << '\n';
		std::cout << "status=EXACT_ALL_FOUR_CATEGORIES\n";
		return 0;
	}
	catch (const std::exception& error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
}
