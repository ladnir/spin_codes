#pragma once

#include <libOTe/Tools/RiffleCode/StructuredSpinB256T128S19Parameters.h>

#include <cryptoTools/Common/Defines.h>

#include <cstddef>

namespace osuCrypto
{
	struct StructuredSpinTestAccess;

	class StructuredSpinB256T128S19
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

			friend class StructuredSpinB256T128S19;
		};

		static constexpr u64 messageBlocks =
			StructuredSpinB256T128S19Parameters::messageBlocks;
		static constexpr u64 codeBlocks =
			StructuredSpinB256T128S19Parameters::codeBlocks;
		static constexpr u64 outerDimension =
			StructuredSpinB256T128S19Parameters::outerDimension;
		static constexpr u64 outerLength =
			StructuredSpinB256T128S19Parameters::outerLength;
		static constexpr u64 outerBlocks = messageBlocks / outerDimension;
		static constexpr u64 epochs =
			StructuredSpinB256T128S19Parameters::epochs;

		StructuredSpinB256T128S19();
		~StructuredSpinB256T128S19();

		StructuredSpinB256T128S19(const StructuredSpinB256T128S19&) = delete;
		StructuredSpinB256T128S19& operator=(const StructuredSpinB256T128S19&) = delete;
		StructuredSpinB256T128S19(StructuredSpinB256T128S19&&) = delete;
		StructuredSpinB256T128S19& operator=(StructuredSpinB256T128S19&&) = delete;

		void init(
			u64 permutationSeed,
			u64 setupCoefficientSeed,
			u64 innerCoefficientSeed,
			bool retainOracleSchedules = false);

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
		static constexpr std::size_t storageBytes = 320;
		alignas(64) std::byte mStorage[storageBytes];

		friend struct StructuredSpinTestAccess;
	};
}
