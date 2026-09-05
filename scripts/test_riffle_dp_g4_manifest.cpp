#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

#include "libOTe/Tools/RiffleCode/RiffleDoubleParityG4.h"

using namespace osuCrypto;

namespace
{
	block basisBlock(u64 bit)
	{
		return bit < 64 ? block(u64{ 1 } << bit, 0) :
			block(0, u64{ 1 } << (bit - 64));
	}

	u64 multiplyX(u64 value)
	{
		const auto high = value >> 63;
		return (value << 1) ^
			(RiffleDoubleParityG4::fieldReduction & (u64{ 0 } - high));
	}

	void applyCoefficientTranspose(
		u64 coefficient,
		const block* __restrict input,
		block* __restrict output)
	{
		for (u64 column = 0; column < 64; ++column)
		{
			block value{};
			auto rows = coefficient;
			while (rows)
			{
				const auto row = static_cast<u64>(std::countr_zero(rows));
				value ^= input[row];
				rows &= rows - 1;
			}
			output[column] = value;
			coefficient = multiplyX(coefficient);
		}
	}

	void referenceInner(
		span<block> word,
		span<const u32> packetOrder)
	{
		constexpr u64 packetWidth = RiffleDoubleParityG4::packetWidth;
		constexpr u64 packetsPerNode = 64 / packetWidth;
		std::array<block, 64> state{};
		std::array<block, 64> result{};
		std::array<block*, 1> resultPointer{ result.data() };
		std::array<const block*, 1> statePointer{ state.data() };

		for (u64 node = RiffleDoubleParityG4::innerNodeCount; node-- > 0;)
		{
			const auto packetBase = node * packetsPerNode;
			for (u64 packet = 0; packet < packetsPerNode; ++packet)
			{
				const auto physical = packetOrder[packetBase + packet];
				std::copy_n(
					word.data() + u64{ physical } * packetWidth,
					packetWidth,
					result.data() + packet * packetWidth);
			}
			statePointer[0] = state.data();
			detail::extendedBch128x64SystematicPacketTransposeStep<64>(
				resultPointer.data(), statePointer.data());
			state = result;
			for (u64 packet = 0; packet < packetsPerNode; ++packet)
			{
				const auto physical = packetOrder[packetBase + packet];
				std::copy_n(
					result.data() + packet * packetWidth,
					packetWidth,
					word.data() + u64{ physical } * packetWidth);
			}
		}
	}

	void referenceOuter(span<const block> word, span<block> output)
	{
		std::array<block, 64> parity0{};
		std::array<block, 64> parity1{};
		std::array<block, 64> message{};
		std::array<block, 64> weighted{};
		ExtendedBch128x64::transpose(
			word.data() + RiffleDoubleParityG4::dataBlockCount * 128,
			parity0.data());
		ExtendedBch128x64::transpose(
			word.data() + (RiffleDoubleParityG4::dataBlockCount + 1) * 128,
			parity1.data());

		u64 coefficient = 1;
		for (u64 outerBlock = 0;
			outerBlock < RiffleDoubleParityG4::dataBlockCount;
			++outerBlock)
		{
			ExtendedBch128x64::transpose(
				word.data() + outerBlock * 128, message.data());
			applyCoefficientTranspose(
				coefficient, parity1.data(), weighted.data());
			for (u64 bit = 0; bit < 64; ++bit)
				output[outerBlock * 64 + bit] =
					message[bit] ^ parity0[bit] ^ weighted[bit];
			coefficient = multiplyX(coefficient);
		}
	}

	bool exactLocalBchKernelCheck()
	{
		std::array<block, 128> word{};
		std::array<block, 64> direct{};
		std::array<block, 64> optimized{};
		for (u64 coordinate = 0; coordinate < 128; ++coordinate)
			word[coordinate] = basisBlock(coordinate);
		ExtendedBch128x64::transpose(word.data(), direct.data());
		ExtendedBch128x64::transposeBlock(word.data(), optimized.data());
		if (direct != optimized)
			return false;

		std::array<block, 64> addend{};
		for (u64 bit = 0; bit < 64; ++bit)
			addend[bit] = basisBlock(bit + 64);
		for (u64 bit = 0; bit < 64; ++bit)
			direct[bit] ^= addend[bit];
		ExtendedBch128x64::transposeAddBlock(
			word.data(), addend.data(), optimized.data());
		return direct == optimized;
	}

	bool exactInnerKernelCheck()
	{
		constexpr u64 packetWidth = RiffleDoubleParityG4::packetWidth;
		constexpr u64 packetCount = 64 / packetWidth;
		std::array<block, 64> hotCurrent{};
		std::array<block, 64> hotState{};
		for (u64 bit = 0; bit < 64; ++bit)
		{
			hotCurrent[bit] = basisBlock(bit);
			hotState[bit] = basisBlock(bit + 64);
		}
		auto directCurrent = hotCurrent;
		auto directState = hotState;
		std::array<block*, packetCount> hotPointers{};
		for (u64 packet = 0; packet < packetCount; ++packet)
			hotPointers[packet] = hotCurrent.data() + packet * packetWidth;
		detail::extendedBch128x64SystematicPacketTransposeHotStateStep<packetWidth>(
			hotPointers.data(), hotState.data());

		std::array<block*, 1> directCurrentPointer{ directCurrent.data() };
		std::array<const block*, 1> directStatePointer{ directState.data() };
		detail::extendedBch128x64SystematicPacketTransposeStep<64>(
			directCurrentPointer.data(), directStatePointer.data());
		return hotCurrent == directCurrent && hotState == directCurrent;
	}

	bool exactCoefficientRecurrenceCheck()
	{
		std::array<block, 64> initial{};
		std::array<block, 64> state{};
		std::array<block, 64> direct{};
		for (u64 bit = 0; bit < 64; ++bit)
			initial[bit] = state[bit] = basisBlock(bit);
		u64 coefficient = 1;
		u64 head = 0;
		for (u64 outerBlock = 0;
			outerBlock < RiffleDoubleParityG4::dataBlockCount;
			++outerBlock)
		{
			applyCoefficientTranspose(
				coefficient, initial.data(), direct.data());
			for (u64 bit = 0; bit < 64; ++bit)
			{
				if (state[(head + bit) & 63] != direct[bit])
					return false;
			}
			const auto tail =
				state[head] ^
				state[(head + 1) & 63] ^
				state[(head + 3) & 63] ^
				state[(head + 4) & 63];
			state[head] = tail;
			head = (head + 1) & 63;
			coefficient = multiplyX(coefficient);
		}
		return true;
	}
}

int main()
{
	if (!exactLocalBchKernelCheck())
	{
		std::cerr << "local BCH transpose kernel mismatch\n";
		return 1;
	}
	if (!exactInnerKernelCheck())
	{
		std::cerr << "packet-four inner kernel mismatch\n";
		return 2;
	}
	if (!exactCoefficientRecurrenceCheck())
	{
		std::cerr << "double-parity coefficient recurrence mismatch\n";
		return 3;
	}

	std::vector<u32> packetOrder(RiffleDoubleParityG4::packetCount);
	for (u64 logical = 0; logical < packetOrder.size(); ++logical)
		packetOrder[logical] = static_cast<u32>(
			(5 * logical + 17) % packetOrder.size());
	std::vector<bool> seen(packetOrder.size());
	for (const auto packet : packetOrder)
	{
		if (seen[packet])
			return 4;
		seen[packet] = true;
	}

	std::vector<block> input(RiffleDoubleParityG4::codeSize);
	u64 state = 0x9e3779b97f4a7c15ULL;
	for (u64 index = 0; index < input.size(); ++index)
	{
		state ^= state >> 12;
		state ^= state << 25;
		state ^= state >> 27;
		const auto low = state * 0x2545f4914f6cdd1dULL;
		state += index + 0xd1b54a32d192ed03ULL;
		const auto high = state ^ (state << 17);
		input[index] = block(low, high);
	}

	auto expectedWord = input;
	referenceInner(expectedWord, packetOrder);
	std::vector<block> expectedOutput(RiffleDoubleParityG4::messageSize);
	referenceOuter(expectedWord, expectedOutput);

	RiffleDoubleParityG4 code;
	code.initOrder(packetOrder);
	auto actualWord = input;
	std::vector<block> actualOutput(RiffleDoubleParityG4::messageSize);
	code.transpose(actualWord, actualOutput);

	if (actualWord != expectedWord)
	{
		std::cerr << "inner transpose mismatch\n";
		return 5;
	}
	if (actualOutput != expectedOutput)
	{
		std::cerr << "outer transpose mismatch\n";
		return 6;
	}
	std::cout << "status=EXACT_IMPLEMENTATION_EQUIVALENCE_VERIFIED\n";
	return 0;
}
