#include <libOTe/Tools/RiffleCode/GolayBa3Rm2SubB240K20.h>
#include <libOTe/Tools/RiffleCode/StructuredSpinB256T128S19.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__linux__)
#include <pthread.h>
#include <sched.h>
#include <unistd.h>
#endif

namespace
{
	using osuCrypto::block;
	using osuCrypto::GolayBa3OuterMode;
	using osuCrypto::GolayBa3Rm2SubB240K20;
	using osuCrypto::StructuredSpinB256T128S19;
	using osuCrypto::u64;

	volatile u64 benchmarkSink = 0;

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

	u64 lowWord(const block& value) noexcept
	{
		return static_cast<u64>(_mm_cvtsi128_si64(value.mData));
	}

	void pinThread(u64 cpu)
	{
#if defined(__linux__)
		cpu_set_t set;
		CPU_ZERO(&set);
		CPU_SET(static_cast<int>(cpu), &set);
		if (pthread_setaffinity_np(pthread_self(), sizeof(set), &set) != 0)
			throw std::runtime_error("failed to pin benchmark thread");
#else
		(void)cpu;
#endif
	}

	void ensureNoOtherBenchmark()
	{
#if defined(__linux__)
		const auto self = static_cast<unsigned long>(getpid());
		for (const auto& entry : std::filesystem::directory_iterator("/proc"))
		{
			const std::string name = entry.path().filename().string();
			if (name.empty() || !std::all_of(name.begin(), name.end(), ::isdigit))
				continue;
			if (std::strtoul(name.c_str(), nullptr, 10) == self)
				continue;
			std::error_code error;
			const auto executable = std::filesystem::read_symlink(entry.path() / "exe", error);
			if (error)
				continue;
			std::string base = executable.filename().string();
			std::transform(base.begin(), base.end(), base.begin(),
				[](unsigned char value) { return static_cast<char>(std::tolower(value)); });
			if (base.find("bench") != std::string::npos)
				throw std::runtime_error("another benchmark executable is active: " + base);
		}
#endif
	}

	double median(std::vector<double> values)
	{
		std::sort(values.begin(), values.end());
		return values[values.size() / 2];
	}

	template<typename Function>
	double timeOnce(Function&& function)
	{
		const auto begin = std::chrono::steady_clock::now();
		function();
		const auto end = std::chrono::steady_clock::now();
		return std::chrono::duration<double, std::milli>(end - begin).count();
	}
}

int main(int argc, char** argv)
{
	try
	{
		ensureNoOtherBenchmark();
		const u64 trials = argc > 1 ? std::strtoull(argv[1], nullptr, 0) : 21;
		if (trials == 0 || (trials & 1) == 0)
			throw std::invalid_argument("trial count must be positive and odd");
		const char* cpuText = std::getenv("SPIN_BENCH_CPU");
		const u64 cpu = cpuText ? std::strtoull(cpuText, nullptr, 0) : 15;
		pinThread(cpu);

		StructuredSpinB256T128S19 baseline;
		baseline.init(
			0x5045524d55544531ULL,
			0x5345545550434f45ULL,
			0x494e4e4552434f45ULL);
		GolayBa3Rm2SubB240K20 reused;
		reused.init(
			0x524f5554452d3234ULL,
			0x42412d332d423234ULL,
			0x524d325355423139ULL,
			GolayBa3OuterMode::Reused);
		GolayBa3Rm2SubB240K20 independent;
		independent.init(
			0x524f5554452d3234ULL,
			0x42412d332d423234ULL,
			0x524d325355423139ULL,
			GolayBa3OuterMode::Independent);

		u64 randomState = 0x42454e43484b3230ULL;
		std::vector<block> baselineInput(StructuredSpinB256T128S19::codeBlocks);
		std::vector<block> finiteInput(GolayBa3Rm2SubB240K20::codeBlocks);
		for (auto& value : baselineInput)
			value = randomBlock(randomState);
		for (auto& value : finiteInput)
			value = randomBlock(randomState);
		std::vector<block> baselineOutput(StructuredSpinB256T128S19::messageBlocks);
		std::vector<block> reusedOutput(GolayBa3Rm2SubB240K20::messageBlocks);
		std::vector<block> independentOutput(GolayBa3Rm2SubB240K20::messageBlocks);
		StructuredSpinB256T128S19::Workspace baselineWorkspace;
		GolayBa3Rm2SubB240K20::Workspace reusedWorkspace;
		GolayBa3Rm2SubB240K20::Workspace independentWorkspace;

		baseline.dualEncodeUnchecked(
			baselineInput.data(), baselineOutput.data(), baselineWorkspace);
		reused.dualEncodeUnchecked(
			finiteInput.data(), reusedOutput.data(), reusedWorkspace);
		independent.dualEncodeUnchecked(
			finiteInput.data(), independentOutput.data(), independentWorkspace);

		const std::string mode = argc > 2 ? argv[2] : "rotated";
		if (mode != "rotated")
		{
			std::vector<double> times;
			times.reserve(trials);
			if (mode == "baseline")
			{
				for (unsigned warm = 0; warm < 3; ++warm)
					baseline.dualEncodeUnchecked(
						baselineInput.data(), baselineOutput.data(), baselineWorkspace);
				for (u64 trial = 0; trial < trials; ++trial)
					times.push_back(timeOnce([&] {
						baseline.dualEncodeUnchecked(
							baselineInput.data(), baselineOutput.data(), baselineWorkspace);
					}));
				benchmarkSink ^= lowWord(baselineOutput[0]);
			}
			else if (mode == "reused")
			{
				for (unsigned warm = 0; warm < 3; ++warm)
					reused.dualEncodeUnchecked(
						finiteInput.data(), reusedOutput.data(), reusedWorkspace);
				for (u64 trial = 0; trial < trials; ++trial)
					times.push_back(timeOnce([&] {
						reused.dualEncodeUnchecked(
							finiteInput.data(), reusedOutput.data(), reusedWorkspace);
					}));
				benchmarkSink ^= lowWord(reusedOutput[0]);
			}
			else if (mode == "independent")
			{
				for (unsigned warm = 0; warm < 3; ++warm)
					independent.dualEncodeUnchecked(
						finiteInput.data(), independentOutput.data(), independentWorkspace);
				for (u64 trial = 0; trial < trials; ++trial)
					times.push_back(timeOnce([&] {
						independent.dualEncodeUnchecked(
							finiteInput.data(), independentOutput.data(), independentWorkspace);
					}));
				benchmarkSink ^= lowWord(independentOutput[0]);
			}
			else
				throw std::invalid_argument("mode must be rotated, baseline, reused, or independent");
			std::cout << std::fixed << std::setprecision(6)
				<< "mode=" << mode
				<< " trials=" << trials
				<< " benchmark_cpu=" << cpu
				<< " median_ms=" << median(times)
				<< " benchmark_sink=0x" << std::hex << benchmarkSink << std::dec
				<< '\n';
			return 0;
		}

		std::vector<double> baselineTimes;
		std::vector<double> reusedTimes;
		std::vector<double> independentTimes;
		baselineTimes.reserve(trials);
		reusedTimes.reserve(trials);
		independentTimes.reserve(trials);
		for (u64 trial = 0; trial < trials; ++trial)
		{
			for (u64 slot = 0; slot < 3; ++slot)
			{
				const u64 which = (trial + slot) % 3;
				if (which == 0)
					baselineTimes.push_back(timeOnce([&] {
						baseline.dualEncodeUnchecked(
							baselineInput.data(), baselineOutput.data(), baselineWorkspace);
					}));
				else if (which == 1)
					reusedTimes.push_back(timeOnce([&] {
						reused.dualEncodeUnchecked(
							finiteInput.data(), reusedOutput.data(), reusedWorkspace);
					}));
				else
					independentTimes.push_back(timeOnce([&] {
						independent.dualEncodeUnchecked(
							finiteInput.data(), independentOutput.data(), independentWorkspace);
					}));
			}
			benchmarkSink ^= lowWord(baselineOutput[trial % baselineOutput.size()]);
			benchmarkSink ^= lowWord(reusedOutput[trial % reusedOutput.size()]);
			benchmarkSink ^= lowWord(independentOutput[trial % independentOutput.size()]);
		}

		const double baselineMedian = median(baselineTimes);
		const double reusedMedian = median(reusedTimes);
		const double independentMedian = median(independentTimes);
		std::cout << std::fixed << std::setprecision(6)
			<< "trials=" << trials << " benchmark_cpu=" << cpu << '\n'
			<< "structured_b256_ms=" << baselineMedian << '\n'
			<< "golay_ba3_b240_reused_ms=" << reusedMedian << '\n'
			<< "golay_ba3_b240_independent_ms=" << independentMedian << '\n'
			<< "reused_vs_structured_percent="
			<< 100.0 * (reusedMedian / baselineMedian - 1.0) << '\n'
			<< "independent_vs_structured_percent="
			<< 100.0 * (independentMedian / baselineMedian - 1.0) << '\n'
			<< "independent_vs_reused_percent="
			<< 100.0 * (independentMedian / reusedMedian - 1.0) << '\n'
			<< "structured_code_blocks=" << StructuredSpinB256T128S19::codeBlocks << '\n'
			<< "finite_code_blocks=" << GolayBa3Rm2SubB240K20::codeBlocks << '\n'
			<< "reused_setup_bytes=" << reused.setupBytes() << '\n'
			<< "independent_setup_bytes=" << independent.setupBytes() << '\n'
			<< "workspace_bytes=" << reused.workspaceBytes() << '\n'
			<< "benchmark_sink=0x" << std::hex << benchmarkSink << std::dec << '\n';
		return 0;
	}
	catch (const std::exception& error)
	{
		std::cerr << "benchmark=ABORT error=" << error.what() << '\n';
		return 1;
	}
}
