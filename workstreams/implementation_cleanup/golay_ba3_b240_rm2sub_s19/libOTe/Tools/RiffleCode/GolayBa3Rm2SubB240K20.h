#pragma once

#include <cryptoTools/Common/Defines.h>

#include <cstddef>

namespace osuCrypto
{
	enum class GolayBa3OuterMode : u8
	{
		Reused = 0,
		Independent = 1,
	};

	struct GolayBa3Rm2SubTestAccess;

	// Transposed encoder for the finite k=2^20, B=240 design point.
	//
	// The parent is a [2,119,680, 1,059,840] binary linear code. The public
	// map returns the first 2^20 parent-message coordinates, which is the
	// transpose corresponding to fixing the final 11,264 message coordinates
	// to zero.
	class GolayBa3Rm2SubB240K20
	{
	public:
		class Workspace
		{
		public:
			Workspace();
			~Workspace();

			Workspace(const Workspace&) = delete;
			Workspace& operator=(const Workspace&) = delete;
			Workspace(Workspace&&) = delete;
			Workspace& operator=(Workspace&&) = delete;

		private:
			static constexpr std::size_t storageBytes = 64;
			alignas(64) std::byte mStorage[storageBytes];

			friend class GolayBa3Rm2SubB240K20;
			friend struct GolayBa3Rm2SubTestAccess;
		};

		static constexpr u64 messageBlocks = u64{ 1 } << 20;
		static constexpr u64 outerLength = 240;
		static constexpr u64 outerDimension = 120;
		static constexpr u64 outerRows = 8'832;
		static constexpr u64 parentMessageBlocks = outerRows * outerDimension;
		static constexpr u64 codeBlocks = outerRows * outerLength;
		static constexpr u64 shortenedMessageBlocks =
			parentMessageBlocks - messageBlocks;
		static constexpr u64 innerStepBlocks = 128;
		static constexpr u64 innerEpochs = codeBlocks / innerStepBlocks;

		static_assert(codeBlocks == 2'119'680);
		static_assert(parentMessageBlocks == 1'059'840);
		static_assert(shortenedMessageBlocks == 11'264);
		static_assert(codeBlocks % innerStepBlocks == 0);

		GolayBa3Rm2SubB240K20();
		~GolayBa3Rm2SubB240K20();

		GolayBa3Rm2SubB240K20(const GolayBa3Rm2SubB240K20&) = delete;
		GolayBa3Rm2SubB240K20& operator=(const GolayBa3Rm2SubB240K20&) = delete;
		GolayBa3Rm2SubB240K20(GolayBa3Rm2SubB240K20&&) = delete;
		GolayBa3Rm2SubB240K20& operator=(GolayBa3Rm2SubB240K20&&) = delete;

		// Seed-derived schedules instantiate the sampled online map. This method
		// does not test the finite note's global G_240 spectrum event; no efficient
		// verifier for that event is supplied by the theory workstream.
		void init(
			u64 routeSeed,
			u64 baSeed,
			u64 innerCoefficientSeed,
			GolayBa3OuterMode outerMode);

		GolayBa3OuterMode outerMode() const;
		u64 setupBytes() const noexcept;
		u64 workspaceBytes() const noexcept;

		void dualEncodeTo(
			const block* __restrict input,
			u64 inputSize,
			block* __restrict output,
			u64 outputSize,
			Workspace& workspace) const;

		void dualEncodeUnchecked(
			const block* __restrict input,
			block* __restrict output,
			Workspace& workspace) const noexcept;

	private:
		static constexpr std::size_t storageBytes = 256;
		alignas(64) std::byte mStorage[storageBytes];

		friend struct GolayBa3Rm2SubTestAccess;
	};
}
