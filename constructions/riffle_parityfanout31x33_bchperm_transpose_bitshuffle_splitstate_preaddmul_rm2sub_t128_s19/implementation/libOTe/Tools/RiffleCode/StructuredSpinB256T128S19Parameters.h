#pragma once

#include <cstdint>

namespace osuCrypto
{
	struct StructuredSpinB256T128S19Parameters
	{
		static constexpr std::uint64_t messageBlocks = std::uint64_t{ 1 } << 20;
		static constexpr std::uint64_t codeBlocks = std::uint64_t{ 1 } << 21;
		static constexpr std::uint64_t blockBits = 128;
		static constexpr std::uint64_t outerDimension = 128;
		static constexpr std::uint64_t outerLength = 256;
		static constexpr std::uint64_t innerStepBlocks = 128;
		static constexpr std::uint64_t innerStateBlocks = 19;
		static constexpr std::uint64_t epochs = messageBlocks / outerDimension;
		static constexpr std::uint64_t innerEpochs = codeBlocks / innerStepBlocks;
		static constexpr std::uint64_t fanoutTargets = 33;
		static constexpr std::uint64_t fanoutSources = 31;
		static constexpr std::uint64_t tileOuterBlocks = 2048;
	};
}
