#include "StructuredSpin_Tests.h"

#include "StructuredSpin_TestSupport.h"

namespace osuCrypto
{
	void StructuredSpin_correctness_test(const oc::CLP&)
	{
		using namespace structured_spin_test;

		verifyBchTranspose();
		verifyRm2SubS19Transpose();

		const auto fixture = makeFixture();
		StructuredSpinB256T128S19 spin;
		spin.init(
			PermutationSeed,
			CoefficientStreamSeed,
			InnerCoefficientSeed,
			true);
		if (StructuredSpinTestAccess::outer(spin).coefficients() != fixture.coefficients)
			throw std::runtime_error("setup coefficient schedule changed");

		StructuredSpinB256T128S19::Workspace workspace;
		std::vector<block> fusedMessage;
		std::vector<block> innerWord;
		verifyCompleteEncoder(
			fixture, spin, workspace, fusedMessage, innerWord);

		std::vector<block> checkedMessage(
			StructuredSpinB256T128S19::messageBlocks);
		spin.dualEncodeTo(
			fixture.source.data(), fixture.source.size(),
			checkedMessage.data(), checkedMessage.size(), workspace);
		if (!equalBlocks(
			checkedMessage.data(), fusedMessage.data(), checkedMessage.size()))
			throw std::runtime_error(
				"checked interface disagrees with fused oracle");
	}
}
