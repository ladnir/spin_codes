#pragma once

#include <array>
#include <bit>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <libOTe/Tools/RiffleCode/StructuredSpinB256T128S19.h>
#include <libOTe/Tools/RiffleCode/StructuredSpin_TestAccess.hpp>

namespace structured_spin_test
{
	using namespace osuCrypto;

	inline constexpr u64 ExpectedChecksum = 0x95c9d722a9539fefULL;
	// The frozen exploration harness XORed this independent legacy schedule
	// self-test into benchmarkSink before hashing the RM2Sub result. Keep the
	// component explicit so the historical checksum remains reproducible.
	inline constexpr u64 FrozenHarnessSelfTestChecksum = 0x9f3b59bea84e65f5ULL;
	inline constexpr u64 PermutationSeed = 0x5354525543543235ULL;
	inline constexpr u64 CoefficientStreamSeed = 0x4649454c4443484bULL;
	inline constexpr u64 InnerCoefficientSeed = 0x524d32434f454646ULL;

	inline u64 splitmix64(u64& state) noexcept
	{
		u64 value = (state += 0x9e3779b97f4a7c15ULL);
		value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
		value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
		return value ^ (value >> 31);
	}

	inline block randomBlock(u64& state) noexcept
	{
		const u64 low = splitmix64(state);
		const u64 high = splitmix64(state);
		return block(low, high);
	}

	inline bool equalBlocks(const block* left, const block* right, u64 size) noexcept
	{
		return std::memcmp(left, right, size * sizeof(block)) == 0;
	}

	inline u64 checksum(const block* values, u64 size) noexcept
	{
		block sum(0, 0);
		for (u64 index = 0; index < size; ++index)
			sum ^= values[index];
		return sum.get<u64>()[0] ^ sum.get<u64>()[1];
	}

	struct Fixture
	{
		std::vector<u64> coefficients;
		std::vector<block> source;
	};

	inline Fixture makeFixture()
	{
		u64 randomState = CoefficientStreamSeed;
		Fixture fixture;
		fixture.coefficients.resize(StructuredSpinB256T128S19::epochs - 1);
		for (auto& coefficient : fixture.coefficients)
		{
			coefficient = splitmix64(randomState);
			coefficient |= static_cast<u64>(coefficient == 0);
		}
		fixture.source.resize(StructuredSpinB256T128S19::codeBlocks);
		for (auto& value : fixture.source)
			value = randomBlock(randomState);
		return fixture;
	}

	inline void verifyBchTranspose()
	{
		u64 randomState = 0x726966666c652d31ULL;
		for (u64 trial = 0; trial < 4; ++trial)
		{
			std::array<block, 256> word;
			std::array<block, 128> optimized;
			std::array<block, 128> reference;
			for (auto& value : word)
				value = randomBlock(randomState);
			ExtendedBch256x128Eq3::transposeBlock(word.data(), optimized.data());
			ExtendedBch256x128Eq3::transposeReference(word.data(), reference.data());
			if (!equalBlocks(optimized.data(), reference.data(), reference.size()))
				throw std::runtime_error("generated BCH transpose disagrees with dense reference");
		}

		std::array<block, 256> word0;
		std::array<block, 256> word1;
		std::array<block, 128> paired0;
		std::array<block, 128> paired1;
		std::array<block, 128> reference0;
		std::array<block, 128> reference1;
		for (auto& value : word0)
			value = randomBlock(randomState);
		for (auto& value : word1)
			value = randomBlock(randomState);
		ExtendedBch256x128Eq3::transposeBlock2(
			word0.data(), word1.data(), paired0.data(), paired1.data());
		ExtendedBch256x128Eq3::transposeBlock(word0.data(), reference0.data());
		ExtendedBch256x128Eq3::transposeBlock(word1.data(), reference1.data());
		if (!equalBlocks(paired0.data(), reference0.data(), reference0.size()) ||
			!equalBlocks(paired1.data(), reference1.data(), reference1.size()))
			throw std::runtime_error("paired BCH transpose disagrees with scalar circuit");
	}

	inline void verifyRm2SubS19Transpose()
	{
		constexpr u64 testEpochs = 5;
		constexpr u64 testBlocks = 128 * testEpochs;
		u64 randomState = 0x524d325355423136ULL;
		std::vector<block> input(testBlocks);
		for (auto& value : input)
			value = randomBlock(randomState);

		RiffleRm2SubS19Transpose inner;
		inner.init(testEpochs, 0x5452414e53504f53ULL);
		std::vector<block> optimized(testBlocks);
		inner.emitReverse(input.data(), testBlocks,
			[&](u64 index, block value) { optimized[index] = value; });

		std::vector<block> reference(testBlocks);
		std::array<block, 19> state;
		state.fill(block(0, 0));
		for (u64 epoch = testEpochs; epoch-- > 0;)
		{
			const block* node = input.data() + 128 * epoch;
			for (u64 point = 0; point < 128; ++point)
			{
				block addend(0, 0);
				u32 column = detail::Rm2Sub19Columns[point];
				while (column)
				{
					const unsigned bit = std::countr_zero(column);
					addend ^= state[bit];
					column &= column - 1;
				}
				reference[128 * epoch + point] = node[point] ^ addend;
			}
			if (epoch == 0)
				break;

			std::array<block, 19> syndrome;
			syndrome.fill(block(0, 0));
			for (u64 point = 0; point < 128; ++point)
			{
				u32 column = detail::Rm2Sub19Columns[point];
				while (column)
				{
					const unsigned bit = std::countr_zero(column);
					syndrome[bit] ^= node[point];
					column &= column - 1;
				}
			}
			if (epoch + 1 == testEpochs)
				state = syndrome;
			else
			{
				std::array<block, 19> multiplied;
				const auto& masks = inner.schedules()[epoch].rowMasks;
				for (unsigned row = 0; row < 19; ++row)
				{
					block value(0, 0);
					u32 mask = masks[row];
					while (mask)
					{
						const unsigned bit = std::countr_zero(mask);
						value ^= state[bit];
						mask &= mask - 1;
					}
					multiplied[row] = value ^ syndrome[row];
				}
				state = multiplied;
			}
		}

		if (!equalBlocks(reference.data(), optimized.data(), testBlocks))
		{
			for (u64 index = 0; index < testBlocks; ++index)
				if (reference[index] != optimized[index])
					throw std::runtime_error(
						"optimized RM2Sub-S19 transpose disagrees with dense reference at block " +
						std::to_string(index));
			throw std::runtime_error(
				"optimized RM2Sub-S19 transpose disagrees with dense reference");
		}
	}

	inline u64 verifyCompleteEncoder(
		const Fixture& fixture,
		StructuredSpinB256T128S19& spin,
		StructuredSpinB256T128S19::Workspace& workspace,
		std::vector<block>& fusedMessage,
		std::vector<block>& innerWord)
	{
		constexpr u64 codeBlocks = StructuredSpinB256T128S19::codeBlocks;
		constexpr u64 messageBlocks = StructuredSpinB256T128S19::messageBlocks;
		constexpr u64 outerBlocks = StructuredSpinB256T128S19::outerBlocks;
		constexpr u64 outerLength = StructuredSpinB256T128S19::outerLength;
		constexpr u64 outerDimension = StructuredSpinB256T128S19::outerDimension;

		std::vector<block> outerWord(codeBlocks);
		std::vector<block> stagedMessage(messageBlocks);
		innerWord.resize(codeBlocks);
		fusedMessage.resize(messageBlocks);

		StructuredSpinTestAccess::inner(spin).emitReverse(fixture.source.data(), codeBlocks,
			[&](u64 index, block value) { innerWord[index] = value; });
		for (u64 innerIndex = 0; innerIndex < codeBlocks; ++innerIndex)
			outerWord[StructuredSpinTestAccess::outer(spin).innerToOuter()[innerIndex]] = innerWord[innerIndex];

		for (u64 outer = 0; outer < outerBlocks; outer += 2)
		{
			for (u64 lane = 0; lane < 2; ++lane)
			{
				block* word = outerWord.data() + (outer + lane) * outerLength;
				const auto& fanout = StructuredSpinTestAccess::outer(spin).parityFanouts()[outer + lane];
				block pivot = word[fanout.targets[0]];
				for (unsigned index = 1; index < fanout.targets.size(); ++index)
					pivot ^= word[fanout.targets[index]];
				for (const u8 sourceCoordinate : fanout.sources)
					word[sourceCoordinate] ^= pivot;
			}
			ExtendedBch256x128Eq3::transposeBlock2(
				outerWord.data() + outer * outerLength,
				outerWord.data() + (outer + 1) * outerLength,
				stagedMessage.data() + outer * outerDimension,
				stagedMessage.data() + (outer + 1) * outerDimension);
		}

		spin.dualEncodeUnchecked(
			fixture.source.data(), fusedMessage.data(), workspace);
		if (!equalBlocks(stagedMessage.data(), fusedMessage.data(), messageBlocks))
			throw std::runtime_error(
				"fused Structured SPIN encoder disagrees with independent staged encoder");

		const u64 canonicalChecksum =
			checksum(fusedMessage.data(), fusedMessage.size()) ^
			checksum(innerWord.data(), innerWord.size());
		const u64 receiptChecksum = canonicalChecksum ^ FrozenHarnessSelfTestChecksum;
		if (receiptChecksum != ExpectedChecksum)
			throw std::runtime_error(
				"Structured SPIN checksum disagrees with frozen receipt: actual=" +
				std::to_string(receiptChecksum) +
				", expected=" + std::to_string(ExpectedChecksum));
		return receiptChecksum;
	}
}
