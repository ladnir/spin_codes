#include <libOTe_Tests/StructuredSpin_TestSupport.h>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#ifdef _WIN32
#include <Windows.h>
#else
#include <pthread.h>
#include <sched.h>
#endif

namespace
{
	using namespace osuCrypto;
	using namespace structured_spin_test;

	double median(std::vector<double> samples)
	{
		std::sort(samples.begin(), samples.end());
		return samples[samples.size() / 2];
	}

	template<typename Function>
	double milliseconds(Function&& function)
	{
		const auto begin = std::chrono::steady_clock::now();
		function();
		const auto end = std::chrono::steady_clock::now();
		return std::chrono::duration<double, std::milli>(end - begin).count();
	}

	template<typename Function>
	std::vector<double> benchmark(u64 trials, Function&& function)
	{
		std::vector<double> samples;
		samples.reserve(trials);
		for (u64 trial = 0; trial < trials; ++trial)
			samples.push_back(milliseconds(function));
		return samples;
	}

	void printSamples(const char* name, const std::vector<double>& samples)
	{
		std::cout << name << "_median_ms=" << median(samples)
			<< ' ' << name << "_samples_ms=";
		for (const auto sample : samples)
			std::cout << sample << ',';
		std::cout << '\n';
	}

	void pinThread(u64 cpu)
	{
#ifdef _WIN32
		if (!SetThreadAffinityMask(GetCurrentThread(), DWORD_PTR{ 1 } << cpu))
			throw std::runtime_error("failed to pin benchmark thread");
		SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_ABOVE_NORMAL);
#else
		cpu_set_t affinity;
		CPU_ZERO(&affinity);
		CPU_SET(cpu, &affinity);
		if (pthread_setaffinity_np(pthread_self(), sizeof(affinity), &affinity) != 0)
			throw std::runtime_error("failed to pin benchmark thread");
#endif
	}
}

int main(int argc, char** argv)
{
	try
	{
		const u64 trials = argc > 1 ? std::strtoull(argv[1], nullptr, 0) : 21;
		if (trials < 3 || !(trials & 1))
			throw std::invalid_argument("trials must be an odd integer at least three");
		const char* cpuText = std::getenv("SPIN_BENCH_CPU");
		const u64 benchmarkCpu = cpuText ? std::strtoull(cpuText, nullptr, 0) : 15;
		pinThread(benchmarkCpu);

		const auto fixture = makeFixture();
		StructuredSpinB256T128S19 spin;
		spin.init(
			PermutationSeed,
			CoefficientStreamSeed,
			InnerCoefficientSeed,
			true);
		StructuredSpinB256T128S19::Workspace workspace;
		std::vector<block> fusedMessage;
		std::vector<block> innerWord;
		const u64 receiptChecksum = verifyCompleteEncoder(
			fixture, spin, workspace, fusedMessage, innerWord);

		const auto integrated = benchmark(trials, [&] {
			spin.dualEncodeUnchecked(
				fixture.source.data(), fusedMessage.data(), workspace);
		});
		const auto materializedInner = benchmark(trials, [&] {
			StructuredSpinTestAccess::inner(spin).emitReverse(fixture.source.data(), fixture.source.size(),
				[&](u64 index, block value) { innerWord[index] = value; });
		});

		std::cout << std::fixed << std::setprecision(6)
			<< "construction=Structured-SPIN-B256-T128-S19"
			<< " benchmark_cpu=" << benchmarkCpu
			<< " trials=" << trials
			<< " correctness=PASS\n";
		printSamples("integrated_full", integrated);
		printSamples("materialized_inner", materializedInner);
		std::cout << "canonical_output_checksum=0x" << std::hex
			<< (receiptChecksum ^ FrozenHarnessSelfTestChecksum) << '\n';
		std::cout << "checksum=0x" << std::hex << receiptChecksum << std::dec << '\n';
		return 0;
	}
	catch (const std::exception& error)
	{
		std::cerr << "error=" << error.what() << '\n';
		return 1;
	}
}
