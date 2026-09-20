#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{
	using u8 = std::uint8_t;
	using u16 = std::uint16_t;
	using u64 = std::uint64_t;

	constexpr u64 dataSymbols = u64{ 1 } << 14;
	constexpr u64 fieldReduction = 0x1b;
	constexpr u64 localGenerator = 0xf4845518b9582a1fULL;
	constexpr u64 emptyKey = 0;

	struct Hash
	{
		static u64 mix(u64 value)
		{
			value ^= value >> 30;
			value *= 0xbf58476d1ce4e5b9ULL;
			value ^= value >> 27;
			value *= 0x94d049bb133111ebULL;
			return value ^ (value >> 31);
		}

		std::size_t operator()(u64 value) const
		{
			return static_cast<std::size_t>(mix(value));
		}
	};

	struct Word128
	{
		u64 low = 0;
		u64 high = 0;
	};

	struct LowWord
	{
		u64 message = 0;
		Word128 codeword{};
		u8 support = 0;
	};

	class MinimumTable
	{
	public:
		explicit MinimumTable(u64 maximumInsertions)
		{
			u64 capacity = 1024;
			while (capacity * 7 < maximumInsertions * 10)
				capacity <<= 1;
			mKeys.resize(capacity);
			mMinimums.assign(capacity, std::numeric_limits<u8>::max());
			mMask = capacity - 1;
		}

		void insert(u64 key, u8 minimum)
		{
			if (key == emptyKey)
				throw std::invalid_argument("zero is not a valid multiplier");
			auto index = Hash::mix(key) & mMask;
			while (mKeys[index] != emptyKey && mKeys[index] != key)
				index = (index + 1) & mMask;
			if (mKeys[index] == emptyKey)
			{
				mKeys[index] = key;
				++mSize;
			}
			mMinimums[index] = std::min(mMinimums[index], minimum);
		}

		u8 find(u64 key) const
		{
			if (key == emptyKey)
				return std::numeric_limits<u8>::max();
			auto index = Hash::mix(key) & mMask;
			while (mKeys[index] != emptyKey)
			{
				if (mKeys[index] == key)
					return mMinimums[index];
				index = (index + 1) & mMask;
			}
			return std::numeric_limits<u8>::max();
		}

		u64 size() const { return mSize; }
		u64 capacity() const { return mKeys.size(); }

	private:
		std::vector<u64> mKeys;
		std::vector<u8> mMinimums;
		u64 mMask = 0;
		u64 mSize = 0;
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
			throw std::invalid_argument("zero has no inverse");
		return power(value, std::numeric_limits<u64>::max() - 1);
	}

	Word128 encode(u64 message)
	{
		Word128 result{};
		while (message)
		{
			const auto bit = static_cast<u64>(std::countr_zero(message));
			result.low ^= localGenerator << bit;
			result.high ^=
				(bit ? localGenerator >> (64 - bit) : 0) | (u64{ 1 } << 63);
			message &= message - 1;
		}
		return result;
	}

	u8 nibbleSupport(Word128 word)
	{
		constexpr u64 lowBits = 0x1111111111111111ULL;
		auto occupiedLow =
			(word.low | (word.low >> 1) | (word.low >> 2) | (word.low >> 3)) & lowBits;
		auto occupiedHigh =
			(word.high | (word.high >> 1) | (word.high >> 2) | (word.high >> 3)) & lowBits;
		return static_cast<u8>(std::popcount(occupiedLow) + std::popcount(occupiedHigh));
	}

	std::array<std::vector<LowWord>, 15> readLowWords(const std::string& path)
	{
		std::ifstream input(path, std::ios::binary);
		if (!input)
			throw std::runtime_error("could not open low-message file");
		std::array<std::vector<LowWord>, 15> result;
		while (true)
		{
			u8 support = 0;
			u64 message = 0;
			input.read(reinterpret_cast<char*>(&support), 1);
			if (!input)
			{
				if (input.eof())
					break;
				throw std::runtime_error("failed reading support byte");
			}
			input.read(reinterpret_cast<char*>(&message), sizeof(message));
			if (!input)
				throw std::runtime_error("partial low-message record");
			if (support >= result.size())
				continue;
			const auto codeword = encode(message);
			if (nibbleSupport(codeword) != support)
				throw std::runtime_error("low message failed independent re-encoding");
			result[support].push_back({ message, codeword, support });
		}
		for (auto& words : result)
		{
			std::sort(words.begin(), words.end(), [](const LowWord& left, const LowWord& right) {
				return left.message < right.message;
			});
		}
		return result;
	}

	bool termLess(const LowWord& left, const LowWord& right)
	{
		return std::pair(left.support, left.message) <
			std::pair(right.support, right.message);
	}
}

int main(int argc, char** argv)
{
	try
	{
		if (argc != 2)
			throw std::invalid_argument("usage: classify_riffle_dp_g4_support38 FILE");
		auto bySupport = readLowWords(argv[1]);
		if (bySupport[11].size() != 20 || bySupport[12].size() != 1'526 ||
			bySupport[13].size() != 37'014 || bySupport[14].empty())
			throw std::runtime_error("input is not a complete local enumeration through support 14");

		std::vector<u64> powers(dataSymbols);
		powers[0] = 1;
		for (u64 exponent = 1; exponent < dataSymbols; ++exponent)
			powers[exponent] = multiplyX(powers[exponent - 1]);

		const u64 repeatedAttempts =
			bySupport[11].size() *
				(bySupport[11].size() + bySupport[12].size() +
				 bySupport[13].size() + bySupport[14].size()) +
			bySupport[12].size() *
				(bySupport[11].size() + bySupport[12].size());
		MinimumTable repeatedCandidates(repeatedAttempts);
		for (u64 parameterSupport = 11; parameterSupport <= 12; ++parameterSupport)
		{
			const auto maximumThirdSupport = 36 - 2 * parameterSupport;
			for (const auto& parameter : bySupport[parameterSupport])
			{
				const auto parameterInverse = inverse(parameter.message);
				for (u64 thirdSupport = 11; thirdSupport <= maximumThirdSupport; ++thirdSupport)
				{
					for (const auto& third : bySupport[thirdSupport])
					{
						const auto coefficient = multiply(third.message, parameterInverse);
						repeatedCandidates.insert(
							coefficient,
							static_cast<u8>(2 * parameterSupport + thirdSupport));
					}
				}
			}
		}

		std::array<u64, 39> oneNew{};
		for (const auto coefficient : powers)
		{
			const auto minimum = repeatedCandidates.find(coefficient);
			if (minimum <= 36)
				++oneNew[minimum];
		}
		std::array<u64, 39> pairP1New{};
		for (u64 left = 0; left + 1 < dataSymbols; ++left)
		{
			for (u64 right = left + 1; right < dataSymbols; ++right)
			{
				const auto minimum = repeatedCandidates.find(powers[left] ^ powers[right]);
				if (minimum <= 36)
					++pairP1New[minimum];
			}
		}

		std::unordered_map<u64, u64, Hash> inverseCache;
		inverseCache.reserve(1'000'000);
		auto cachedInverse = [&](u64 value) {
			const auto found = inverseCache.find(value);
			if (found != inverseCache.end())
				return found->second;
			const auto result = inverse(value);
			inverseCache.emplace(value, result);
			return result;
		};

		std::unordered_map<u64, u8, Hash> additiveRatioMinimum;
		additiveRatioMinimum.reserve(1'000'000);
		std::array<u64, 39> additiveUnorderedTriples{};
		for (u64 leftSupport = 11; leftSupport <= 12; ++leftSupport)
		{
			for (u64 rightSupport = leftSupport; rightSupport <= 13; ++rightSupport)
			{
				const auto& leftWords = bySupport[leftSupport];
				const auto& rightWords = bySupport[rightSupport];
				for (u64 leftIndex = 0; leftIndex < leftWords.size(); ++leftIndex)
				{
					const auto rightBegin =
						leftSupport == rightSupport ? leftIndex + 1 : 0;
					for (u64 rightIndex = rightBegin; rightIndex < rightWords.size(); ++rightIndex)
					{
						const auto& left = leftWords[leftIndex];
						const auto& right = rightWords[rightIndex];
						LowWord third{
							left.message ^ right.message,
							{ left.codeword.low ^ right.codeword.low,
							  left.codeword.high ^ right.codeword.high },
							0
						};
						third.support = nibbleSupport(third.codeword);
						std::array<LowWord, 3> terms{ left, right, third };
						std::sort(terms.begin(), terms.end(), termLess);
						if (terms[0].message != left.message ||
							terms[1].message != right.message)
							continue;
						const auto total = static_cast<u64>(
							terms[0].support + terms[1].support + terms[2].support);
						if (total < 33 || total > 38)
							continue;
						++additiveUnorderedTriples[total];
						for (u64 numerator = 0; numerator < 3; ++numerator)
						{
							for (u64 denominator = 0; denominator < 3; ++denominator)
							{
								if (numerator == denominator)
									continue;
								const auto ratio = multiply(
									terms[numerator].message,
									cachedInverse(terms[denominator].message));
								auto [position, inserted] =
									additiveRatioMinimum.emplace(ratio, static_cast<u8>(total));
								if (!inserted)
									position->second = std::min(position->second, static_cast<u8>(total));
							}
						}
					}
				}
			}
		}

		std::array<u64, 39> pairP0New{};
		for (u64 gap = 1; gap < dataSymbols; ++gap)
		{
			const auto found = additiveRatioMinimum.find(powers[gap]);
			if (found != additiveRatioMinimum.end())
				pairP0New[found->second] += dataSymbols - gap;
		}

		std::array<u64, 39> tripleNew{};
		for (u64 gap2 = 2; gap2 < dataSymbols; ++gap2)
		{
			const auto denominator = powers[gap2] ^ 1;
			const auto denominatorInverse = inverse(denominator);
			const auto gap2Ratio = multiply(powers[gap2], denominatorInverse);
			auto gap1Ratio = denominatorInverse;
			for (u64 gap1 = 1; gap1 < gap2; ++gap1)
			{
				gap1Ratio = multiplyX(gap1Ratio);
				const auto ratio = gap1Ratio ^ gap2Ratio;
				const auto found = additiveRatioMinimum.find(ratio);
				if (found != additiveRatioMinimum.end())
					tripleNew[found->second] += dataSymbols - gap2;
			}
		}

		std::cout << "candidate=Riffle DP g=4\n";
		std::cout << "local_support_11=" << bySupport[11].size() << '\n';
		std::cout << "local_support_12=" << bySupport[12].size() << '\n';
		std::cout << "local_support_13=" << bySupport[13].size() << '\n';
		std::cout << "local_support_14=" << bySupport[14].size() << '\n';
		std::cout << "repeated_candidate_coefficients=" << repeatedCandidates.size() << '\n';
		std::cout << "repeated_candidate_capacity=" << repeatedCandidates.capacity() << '\n';
		std::cout << "additive_ratio_classes=" << additiveRatioMinimum.size() << '\n';
		for (u64 support = 33; support <= 38; ++support)
		{
			std::cout << "support=" << support;
			if (support <= 36)
				std::cout << " one_new_dirty_classes=" << oneNew[support]
					<< " pair_p1_new_dirty_classes=" << pairP1New[support];
			else
				std::cout << " repeated_categories=NOT_CLASSIFIED";
			std::cout << " pair_p0_new_dirty_classes=" << pairP0New[support]
				<< " triple_new_dirty_classes=" << tripleNew[support]
				<< " additive_unordered_triples=" << additiveUnorderedTriples[support]
				<< '\n';
		}
		std::cout << "one_exact_through=36\n";
		std::cout << "pair_p1_exact_through=36\n";
		std::cout << "pair_p0_exact_through=38\n";
		std::cout << "triple_exact_through=38\n";
		std::cout << "status=EXACT_SUPPORT38_CLASSIFICATION\n";
		return 0;
	}
	catch (const std::exception& error)
	{
		std::cerr << "error: " << error.what() << '\n';
		return 1;
	}
}
