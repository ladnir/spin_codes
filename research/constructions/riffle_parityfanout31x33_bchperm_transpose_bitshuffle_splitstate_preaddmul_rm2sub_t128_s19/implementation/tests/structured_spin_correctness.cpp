#include <libOTe_Tests/StructuredSpin_TestSupport.h>

#include <iomanip>
#include <iostream>

int main()
{
	try
	{
		using namespace osuCrypto;
		using namespace structured_spin_test;

		verifyBchTranspose();
		std::cout << "bch_dense=PASS\n";
		verifyRm2SubS19Transpose();
		std::cout << "rm2sub_s19_dense=PASS\n";

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
		const u64 receiptChecksum = verifyCompleteEncoder(
			fixture, spin, workspace, fusedMessage, innerWord);
		std::vector<block> checkedMessage(StructuredSpinB256T128S19::messageBlocks);
		spin.dualEncodeTo(
			fixture.source.data(), fixture.source.size(),
			checkedMessage.data(), checkedMessage.size(), workspace);
		if (!equalBlocks(
			checkedMessage.data(), fusedMessage.data(), checkedMessage.size()))
			throw std::runtime_error("checked interface disagrees with fused oracle");
		const u64 canonicalChecksum = receiptChecksum ^ FrozenHarnessSelfTestChecksum;

		std::cout << "staged_complete=PASS\n"
			<< "fused_complete=PASS\n"
			<< "checked_interface=PASS\n"
			<< "canonical_output_checksum=0x" << std::hex << canonicalChecksum << '\n'
			<< "frozen_harness_selftest_checksum=0x" << FrozenHarnessSelfTestChecksum << '\n'
			<< "checksum=0x" << std::hex << receiptChecksum << std::dec << '\n'
			<< "correctness=PASS\n";
		return 0;
	}
	catch (const std::exception& error)
	{
		std::cerr << "error=" << error.what() << '\n';
		return 1;
	}
}
