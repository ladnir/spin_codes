#include <libOTe/Tools/RiffleCode/GolayBa3Rm2SubB240K20.h>
#include <libOTe/Tools/RiffleCode/GolayBa3Rm2Sub_TestAccess.hpp>

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{
	using osuCrypto::block;
	using osuCrypto::GolayBa3OuterMode;
	using osuCrypto::GolayBa3Rm2SubB240K20;
	using osuCrypto::GolayBa3Rm2SubTestAccess;
	using osuCrypto::u64;

	u64 splitmix64(u64& state) noexcept
	{
		u64 value = (state += 0x9e3779b97f4a7c15ULL);
		value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
		value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
		return value ^ (value >> 31);
	}

	block randomBlock(u64& state) noexcept
	{
		return block(splitmix64(state), splitmix64(state));
	}

	u64 hashBlocks(const block* values, u64 count) noexcept
	{
		u64 hash = 0xcbf29ce484222325ULL;
		for (u64 index = 0; index < count; ++index)
		{
			alignas(16) u64 words[2];
			_mm_store_si128(reinterpret_cast<__m128i*>(words), values[index].mData);
			for (u64 word : words)
			{
				hash ^= word;
				hash *= 0x100000001b3ULL;
				hash ^= hash >> 32;
			}
		}
		return hash;
	}

	void verifyMode(GolayBa3OuterMode mode, u64 expectedHash)
	{
		GolayBa3Rm2SubB240K20 code;
		code.init(
			0x524f5554452d3234ULL,
			0x42412d332d423234ULL,
			0x524d325355423139ULL,
			mode);
		GolayBa3Rm2SubTestAccess::validateSetup(code);

		u64 randomState = 0x46494e4954454b32ULL;
		std::vector<block> input(GolayBa3Rm2SubB240K20::codeBlocks);
		for (auto& value : input)
			value = randomBlock(randomState);
		std::vector<block> optimized(GolayBa3Rm2SubB240K20::messageBlocks);
		std::vector<block> reference(GolayBa3Rm2SubB240K20::messageBlocks);
		GolayBa3Rm2SubB240K20::Workspace optimizedWorkspace;
		GolayBa3Rm2SubB240K20::Workspace referenceWorkspace;

		code.dualEncodeTo(
			input.data(), input.size(), optimized.data(), optimized.size(),
			optimizedWorkspace);
		GolayBa3Rm2SubTestAccess::dualEncodeReference(
			code, input.data(), reference.data(), referenceWorkspace);
		if (!std::equal(optimized.begin(), optimized.end(), reference.begin()))
			throw std::runtime_error("optimized encoder disagrees with scalar staged reference");

		const u64 hash = hashBlocks(optimized.data(), optimized.size());
		if (expectedHash != 0 && hash != expectedHash)
			throw std::runtime_error("deterministic correctness checksum changed");
		std::cout
			<< "mode=" << (mode == GolayBa3OuterMode::Reused ? "reused" : "independent")
			<< " checksum=0x" << std::hex << std::setw(16) << std::setfill('0') << hash
			<< std::dec
			<< " setup_bytes=" << code.setupBytes()
			<< " workspace_bytes=" << code.workspaceBytes()
			<< '\n';
	}
}

int main()
{
	try
	{
		verifyMode(GolayBa3OuterMode::Reused, 0xc4108a1282911d32ULL);
		verifyMode(GolayBa3OuterMode::Independent, 0x2bd00619556f9a86ULL);
		std::cout << "golay_ba3_b240_k20_correctness=PASS\n";
		return 0;
	}
	catch (const std::exception& error)
	{
		std::cerr << "golay_ba3_b240_k20_correctness=FAIL error=" << error.what() << '\n';
		return 1;
	}
}
