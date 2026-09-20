#include "StructuredSpinB256T128S19.h"

#include "RiffleExactPermFieldCheckpoint.h"
#include "StructuredSpin_TestAccess.hpp"

#include <memory>
#include <new>
#include <stdexcept>

namespace osuCrypto
{
	namespace
	{
		struct StructuredSpinImplementation
		{
			RiffleExactPermFieldCheckpoint outer;
			RiffleRm2SubS19Transpose inner;
			bool initialized = false;
		};

		using StructuredSpinWorkspace =
			RiffleExactPermFieldCheckpoint::Workspace;

		static_assert(sizeof(StructuredSpinImplementation) <= 320);
		static_assert(alignof(StructuredSpinImplementation) <= 64);
		static_assert(sizeof(StructuredSpinWorkspace) <= 64);
		static_assert(alignof(StructuredSpinWorkspace) <= 64);

		static_assert(RiffleExactPermFieldCheckpoint::messageBlocks ==
			StructuredSpinB256T128S19Parameters::messageBlocks);
		static_assert(RiffleExactPermFieldCheckpoint::codeBlocks ==
			StructuredSpinB256T128S19Parameters::codeBlocks);
		static_assert(RiffleExactPermFieldCheckpoint::outerDimension ==
			StructuredSpinB256T128S19Parameters::outerDimension);
		static_assert(RiffleExactPermFieldCheckpoint::outerLength ==
			StructuredSpinB256T128S19Parameters::outerLength);
		static_assert(RiffleExactPermFieldCheckpoint::epochs ==
			StructuredSpinB256T128S19Parameters::epochs);
		static_assert(RiffleRm2SubS19Transpose::stepBlocks ==
			StructuredSpinB256T128S19Parameters::innerStepBlocks);
		static_assert(RiffleExactPermFieldCheckpoint::codeBlocks /
			RiffleRm2SubS19Transpose::stepBlocks ==
			StructuredSpinB256T128S19Parameters::innerEpochs);
		static_assert(RiffleRm2SubS19Transpose::stateBlocks ==
			StructuredSpinB256T128S19Parameters::innerStateBlocks);

		StructuredSpinImplementation& implementation(std::byte* storage) noexcept
		{
			return *std::launder(
				reinterpret_cast<StructuredSpinImplementation*>(storage));
		}

		const StructuredSpinImplementation& implementation(
			const std::byte* storage) noexcept
		{
			return *std::launder(
				reinterpret_cast<const StructuredSpinImplementation*>(storage));
		}

		StructuredSpinWorkspace& workspace(std::byte* storage) noexcept
		{
			return *std::launder(
				reinterpret_cast<StructuredSpinWorkspace*>(storage));
		}
	}

	StructuredSpinB256T128S19::Workspace::Workspace()
	{
		::new (static_cast<void*>(mStorage)) StructuredSpinWorkspace;
	}

	StructuredSpinB256T128S19::Workspace::~Workspace()
	{
		std::destroy_at(&workspace(mStorage));
	}

	StructuredSpinB256T128S19::StructuredSpinB256T128S19()
	{
		::new (static_cast<void*>(mStorage)) StructuredSpinImplementation;
	}

	StructuredSpinB256T128S19::~StructuredSpinB256T128S19()
	{
		std::destroy_at(&implementation(mStorage));
	}

	void StructuredSpinB256T128S19::init(
		u64 permutationSeed,
		u64 setupCoefficientSeed,
		u64 innerCoefficientSeed,
		bool retainOracleSchedules)
	{
		auto& state = implementation(mStorage);
		state.outer.init(
			permutationSeed, setupCoefficientSeed, retainOracleSchedules);
		state.inner.init(
			codeBlocks / RiffleRm2SubS19Transpose::stepBlocks,
			innerCoefficientSeed);
		state.initialized = true;
	}

	void StructuredSpinB256T128S19::dualEncodeTo(
		const block* __restrict input,
		u64 inputSize,
		block* __restrict output,
		u64 outputSize,
		Workspace& encodeWorkspace) const
	{
		if (!implementation(mStorage).initialized)
			throw std::logic_error("StructuredSpinB256T128S19 is not initialized");
		if (inputSize != codeBlocks || outputSize != messageBlocks)
			throw std::invalid_argument("Structured SPIN span size mismatch");
		dualEncodeUnchecked(input, output, encodeWorkspace);
	}

	void StructuredSpinB256T128S19::dualEncodeUnchecked(
		const block* __restrict input,
		block* __restrict output,
		Workspace& encodeWorkspace) const noexcept
	{
		const auto& state = implementation(mStorage);
		state.outer.dualEncodePackedInnerUnchecked<32>(
			input, output, workspace(encodeWorkspace.mStorage), state.inner);
	}

	const RiffleExactPermFieldCheckpoint& StructuredSpinTestAccess::outer(
		const StructuredSpinB256T128S19& spin) noexcept
	{
		return implementation(spin.mStorage).outer;
	}

	const RiffleRm2SubS19Transpose& StructuredSpinTestAccess::inner(
		const StructuredSpinB256T128S19& spin) noexcept
	{
		return implementation(spin.mStorage).inner;
	}
}
