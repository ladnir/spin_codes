// Reproduce exact input/output slice rows B[t,*] for the committed
// extended-BCH [128,64,22] cyclic encoder.  One enumeration simultaneously
// produces t and 64-t by XORing with the all-input word.

#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <vector>

namespace
{
	constexpr std::uint64_t generator = 0xf4845518b9582a1fULL;
	std::array<std::uint64_t, 64> lowRows{};
	std::array<std::uint64_t, 64> highRows{};
	std::uint64_t allLow = 0;
	std::uint64_t allHigh = 0;
	using Histogram = std::array<std::uint64_t, 129>;

	struct Task
	{
		unsigned start;
		unsigned remaining;
		std::uint64_t low;
		std::uint64_t high;
	};

	void enumerate(
		unsigned start,
		unsigned remaining,
		std::uint64_t low,
		std::uint64_t high,
		Histogram& lowHistogram,
		Histogram& highHistogram)
	{
		if (remaining == 0)
		{
			++lowHistogram[std::popcount(low) + std::popcount(high)];
			++highHistogram[
				std::popcount(low ^ allLow) + std::popcount(high ^ allHigh)];
			return;
		}
		const auto last = 64u - remaining;
		for (auto index = start; index <= last; ++index)
				enumerate(
					index + 1,
					remaining - 1,
					low ^ lowRows[index],
					high ^ highRows[index],
					lowHistogram,
					highHistogram);
	}

	void print(unsigned inputWeight, const Histogram& histogram)
	{
		for (unsigned outputWeight = 0; outputWeight <= 128; ++outputWeight)
			if (histogram[outputWeight])
				std::cout << inputWeight << ',' << outputWeight << ','
					<< histogram[outputWeight] << '\n';
	}
}

int main(int argc, char** argv)
{
	const auto weight = argc > 1 ? static_cast<unsigned>(std::strtoul(argv[1], nullptr, 0)) : 6;
	const auto requestedThreads = argc > 2
		? static_cast<unsigned>(std::strtoul(argv[2], nullptr, 0))
		: std::thread::hardware_concurrency();
	if (weight > 32)
	{
		std::cerr << "weight must be at most 32\n";
		return 1;
	}
	for (unsigned row = 0; row < 64; ++row)
	{
		lowRows[row] = generator << row;
		highRows[row] = (row ? generator >> (64 - row) : 0) | (std::uint64_t{1} << 63);
		allLow ^= lowRows[row];
		allHigh ^= highRows[row];
	}

	std::vector<Task> tasks;
	if (weight == 0)
		tasks.push_back({0, 0, 0, 0});
	else if (weight == 1)
		for (unsigned first = 0; first < 64; ++first)
			tasks.push_back({first + 1, 0, lowRows[first], highRows[first]});
	else
		for (unsigned first = 0; first + weight <= 64; ++first)
			for (unsigned second = first + 1; second + weight - 1 <= 64; ++second)
				tasks.push_back({
					second + 1,
					weight - 2,
					lowRows[first] ^ lowRows[second],
					highRows[first] ^ highRows[second]});

	const auto threadCount = std::max(
		1u, std::min(requestedThreads ? requestedThreads : 1u, static_cast<unsigned>(tasks.size())));
	std::vector<Histogram> lowHistograms(threadCount);
	std::vector<Histogram> highHistograms(threadCount);
	std::atomic<std::size_t> nextTask = 0;
	std::vector<std::thread> workers;
	const auto started = std::chrono::steady_clock::now();
	for (unsigned threadIndex = 0; threadIndex < threadCount; ++threadIndex)
		workers.emplace_back([&, threadIndex] {
			while (true)
			{
				const auto taskIndex = nextTask.fetch_add(1, std::memory_order_relaxed);
				if (taskIndex >= tasks.size())
					break;
				const auto& task = tasks[taskIndex];
				enumerate(
					task.start,
					task.remaining,
					task.low,
					task.high,
					lowHistograms[threadIndex],
					highHistograms[threadIndex]);
			}
		});
	for (auto& worker : workers)
		worker.join();

	Histogram lowHistogram{};
	Histogram highHistogram{};
	for (unsigned threadIndex = 0; threadIndex < threadCount; ++threadIndex)
		for (unsigned outputWeight = 0; outputWeight <= 128; ++outputWeight)
		{
			lowHistogram[outputWeight] += lowHistograms[threadIndex][outputWeight];
			highHistogram[outputWeight] += highHistograms[threadIndex][outputWeight];
		}
	const auto elapsed = std::chrono::duration<double>(
		std::chrono::steady_clock::now() - started).count();
	std::cerr << "enumerated C(64," << weight << ") with " << threadCount
		<< " threads in " << elapsed << " seconds\n";
	std::cout << "input_weight,output_weight,count\n";
	print(weight, lowHistogram);
	if (weight != 32)
		print(64 - weight, highHistogram);
}
